import sys, os, shutil

cfg_fn = "run_config.json"

configfile:
    cfg_fn

workdir:
    config["workdir"]

# include software/program information
# HACK: to get current Snakefile directory
_pipeline_dir = os.path.abspath(sys.path[0])
include:
    _pipeline_dir + "/toolsinfo"

o_dir = conifg["workdir"]
d_dir =  config ["data_dir"]
delivery = config["delivery"]

rule all:
    input: "report.html"

rule make_key:
    input: 'run_config.json'
    output: "map_key.txt"
    message: "make the map key file from information in the run_config"
    run:
        f_out = open(wildcards.output, 'w')
        f_out.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimerSequence\tSequenceType\tGroup\tDescription\n")
        for treatment, samples in sorted(config["group"].items()):
            for sampleID in samples:
                seq_type = config.get("sequence_type", config["sample_sequence_type"][sampleID])
                f_out.write("\t".join([sampleID, "", config["F_primer"], config["R_primer"], seq_type, treatment, sampleID]))
                f_out.write("\n")

rule preprocess:
    input: lambda wc: os.path.join(config["data_dir"],config[wildcards.sample])
    output: "seq/{sample}_merged.fastq"
    threads: 8
    message: "Convert unaligned BAM file, {input}, to fasta file"
    run:
        seq_type = config.get("sequence_type", config["sample_sequence_type"][sample])
        if seq_type == "single":
            for fn in input:
                namebits =  os.path.basename(fn)
                if namebits[-1]  == "bam": # PGM unaligned bam file
                    shell("{BAM2FASTX} -a -Q --all {fn} >> {output}")
                elif namebits[-1].lower() == "sff":
                    shell("{QIIME}/process_sff.py -f -i {fn} -o {output}")
        elif seq_type == "paired":
            min_overlap = config.get('paired_min_overlap', 10),
            max_overlap = config.get('paired_max_overlap', 250)
            for i, fn in enumerate([f for f in input if "R1" in f]):
                r1 = fn
                r2 = fn.replace("R1", "R2")
                if not os._exists(r2): print("ERROR: no R2 reads for %s" % (r1))

                # use flash to merge paired rads
                flash_log = "logs/flash_merge_{}.log".format(wildcards.sample )
                flash_out = o_dir + "./amplicon_fastq/{wildcards.sample}_%d/" % ( d+1, )

                # trim illumina sequencing primers:
                shell("{SKEWER} -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
                 "-y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT "
                 "-k 15 -o {flash_out} -t {threads} {input}")

                r1_trim = {flash_out} + "pair1.fastq.gz"
                r2_trim = {flash_out} + "pair2.fastq.gz"

                # shell("{CUTADAPT} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {r1_trim} {r1}")
                # shell("{CUTADAPT} -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o {r2_trim} {r2}")

                shell("{FLASH}  {r1_trim} {r2_trim} -m {min_overlap} -M {max_overlap} -d {flash_out} >> {flash_log} ")
                shell("cat {flash_out}/out.extendedFrags.fastq >> {output}")

def reverse_complement(seq):
    complement = dict(zip("ATGCYRWSKMDVHBN", "TACGRYWSMKHBDVN"))#{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq.upper())
    bases = reversed([complement.get(base, "N") for base in bases])
    bases = ''.join(bases)
    return bases

rule primer_trim:
    input:  "seq/{sample}_merged.fastq"
    output: "seq/{sample}_trimmed.fastq.gz"
    params:
        log = "logs/cut_primer_{wildcards.sample}.log",
        for = lambda config["F_primer"],
        rev = reverse_complement(config["R_primer"]),
        min_length = config["min_amplicon_length"] - len(for) - len(rev)
    threads: 8
    message: "  Trimming amplicon primer sequences..."
    shell: '{CUTADAPT} -g Forward={params.for} -a Reverse={params.rev} -k 15 -m {min_length} -o {output} {input} 2>&1 > {params.log}'

rule split_libraries:
    input: ["amplicon_fastq/{sample}.fastq" for samples in config["treatments"] for sample in samples]
    output: "sl_libraries/seqs.fna"
    params:
        sample_fq_list=",".join(wildcards.input)
        sample_id_list=",".join([sample for samples in config["treatments"] for sample in samples])
        o_dir = "sl_libraries"
    message: "Running split_libraries..."
    shell: 'head -2 map_key.txt > {o_dir}/one_key.txt;' \
           'source /mnt/software/qiime/activate.sh;'
           'split_libraries.py -i {params.sample_fq_list} --sample_id {params.sample_id_list} -o {params.fna_out}'
           '-m {o_dir}/one_key.txt -q 19 -l {min_amplicon_length} -L {max_amplicon_length} --barcode_type "not-barcoded";' \
           'cp {params.o_dir}/split_library_log.txt logs'


rule OTU_picking:
    input: "sl_libraries/seqs_primer_cut.fna"
    output: "otus/otus_sorted.biom", "otus/rep_set.tre", "", "", ""
    params:
        ref  = config.get("reference", "denovo")
        par  = "otus/otu_picking_params.txt"
        o_dir = "/otus",
        hom  = config.get("homology", "97")
        its  = {ITS},
        r16S = {GreenGene},
        r18S = {RDP}
    log: "logs/otu_picking.log"
    threads: 8
    message: "OTU Picking"
    run:
        if ({params.ref} == "denovo"):
            shell("pick_de_novo_otus.py -i {input} -o {params.otus} -f -a -O {threads}")
        elif ({params.ref} == "ITS"):
            opt_params="assign_taxonomy:id_to_taxonomy_fp {ITS}/taxonomy/{params.hom}_otu_taxonomy.txt\n"
                       "assign_taxonomy:reference_seqs_fp {ITS}/rep_sep/{params.hom}_otus.fasta\n""
            shell("cat {_pipeline_dir}/params_base.txt > {params.par}")
            with open(params.par, 'wa') as f_out:
                f_out.write("assign_taxonomy:id_to_taxonomy_fp {ITS}/taxonomy/{params.hom}_otu_taxonomy.txt\n"
                  "assign_taxonomy:reference_seqs_fp {ITS}/rep_sep/{params.hom}_otus.fasta\n")
            shell("pick_open_reference_otus.py -i {input} -o {output.otus} -r {params.ref} "
                  "-p {params.par} -f -a -O {threads}")
        elif ({params.ref} == "16S"):
            shell("pick_open_reference_otus.py -i {input} -o {output.otus} -r {params.ref} -f -a -O {threads}")
        elif ({params.ref} == "18S"):
            pass
        else:
            shell("source /mnt/software/qiime/activate.sh;"
                  "pick_open_reference_otus.py -i {input} -o {output.o_dir} -r {params.ref} -f -a -O {threads};"
                  "cat {output.o_dir}/log_*.txt > {log}")
        shell("source /mnt/software/qiime/activate.sh;"
              "sort_otu_table.py -i {output.o_dir}/otu_table_mc2_w_tax.biom -s SampleID -m {map} -o {output[0]}")


rule biom_summary:
    input: "otus/otus_sorted.biom"
    output: "otus/biom_summary.txt"
    message: "Creating biom_table_summary..."
    shell: "biom summarize-table -i {input} -o {output}"

rule OTU_heatmap:
    input:
        biom = "otus/otus_sorted.biom",
        map = "map_key.txt"
    output: "heatmap/otus_sorted.html"
    params:
        o_dir = "heatmap"
    message: "Making OTU heatmap..."
    shell: "source /mnt/software/qiime/activate.sh;"
           "make_otu_heatmap_html.py -i {input.biom} -m {input.map} -o {params.o_dir} -n 2"

rule OTU_network:
    input:
        map = "map_key.txt",
        biom = "otus/otus_sorted.biom"
    output: "network/real_edge_table.txt"
    params:
        o_dir = "."
    message: "Creating OTU network..."
    shell: "source /mnt/software/qiime/activate.sh;"
           "make_otu_network.py -m {input.map} -i {input.biom} -o {params.o_dir}"

rule sum_taxa:
    input:
        map = "map_key.txt",
        biom = "otus/otus_sorted.biom"
    output:
        "taxa_summary/taxa_summary_plots/bar_charts.html"
    params:
        o_dir = "taxa_summary"
    log: "logs/taxa_summary.log"
    message: "Creating taxa summary plots..."
    shell: "source /mnt/software/qiime/activate.sh;"
           "summarize_taxa_through_plots.py -f -i {input.biom} -m {input.map} -o {params.o_dir} "
            "-p /mnt/grl/brc/application/qiime_pipeline.0.1/summarize_taxa_params.txt;" \
            "cat taxa_summary/log_*.txt > {log}"

rule alpha_rarefaction:
    input:
        map = "map_key.txt",
        biom = "otus/otus_sorted.biom"
        tre  = "otus/rep_set.tre"
    output:  'arare/alpha_rarefaction_plots/rarefaction_plots.html'
    params:
        o_dir = "arare"
    log: "logs/alpha_rarefaction.log"
    threads: 8
    message: "Creating alpha rarefaction plots..."
    shell: "source /mnt/software/qiime/activate.sh;"
           "alpha_rarefaction.py -f -i {input.biom} -m {input.map} "
            "-o {params.o_dir} -p /mnt/grl/brc/application/qiime_pipeline.0.1/alpha_params.txt -t "
            "{input.tre'} -a -O {threads};" \
            "cat {params.o_dir}/log_*.txt > {log}"


rule beta_diversity:
    # use the minimum reads number as the parameter e for even sampling
    input:
        map = "map_key.txt",
        biom = "otus/otus_sorted.biom"
        tree = "otus/rep_set.tre"
        biom_summary = "otus/biom_summary.txt"
    output: "bdiv/bray_curtis_emperor_pcoa_plot/index.html"
    params:
        o_dir = "bdiv"
    log: "logs/beta_diversity.log"
    threads: 8
    message: "Creating beta diversity plots..."
    shell: 'line=`grep " Min:" {biom_summary}`; min=${{line:6}}; min=(${{min//./ }});' \
           "source /mnt/software/qiime/activate.sh;"
           'beta_diversity_through_plots.py -f -i {input.biom} -m {input.map} -t {input.tree} ' \
           '-e $min -a -O {threads} -o {params.o_dir};' \
           'cat {params.o_dir}/log_*.txt > {log}'

rule report:
    input: "bdiv/bray_curtis_emperor_pcoa_plot/index.html", \
           'arare/alpha_rarefaction_plots/rarefaction_plots.html', \
           "taxa_summary/taxa_summary_plots/bar_charts.html", \
           "network/real_edge_table.txt",\
           "heatmap/otus_sorted.html",\
           "otus/biom_summary.txt",\
           "otus/otus_sorted.biom",\
           "otus/rep_set.tre"
    output: "report.html"
    shell: 'echo "place holder" > report.html'