QIIME pipeline
==============

Clone the pipeline

  git clone file:///mnt/grl/brc/version_control/repositories/qiime_sm.git/

log in (CentOS 6 or above) with X11 enabled (use NX client or ssh -Y)
(QIIME only works in CentOS 6 with X windows)

Without X11, the bug in QIIME plotting functions makes it crush
Make sure /mnt/software/anaconda/bin/ is in your PATH so the proper python environment 
can be invoked:

  # Activate python virtual environment qiime1.9
  source activate qiime1.9

Create a run config file according to the run_config.json in the pipeline.

Run QIIME pipeline script, subsititute path and the run_cfg with yours:

  path/to/qiime_sm/run_qiime.sh run_cfg.json

The first few step can be run on SGE to take advantage of the cluster (reads trimming, 
join, OTU assignment etc.).

  # otu_picking instruct the pipeline to run until OTUs are assigned
  # -j specify how many SGE nodes can be used. 
  # -p prints out the command whenever appropriate
  path/to/qiime_sm/run_qiime_sge.sh run_cfg.json otu_picking -j 5 -p

Then finish the rest with run_qiime.sh locally under X windows environment.

To get help with options:

  path/to/qiime_sm/run_qiime.sh -h

Run configuration file
======================

A run config file in JSON format specifies the data path, work directory, reference, 
and other important parameters for the QIIME data analysis.

An example is:

run_config.json

 {
   "PIname": "Bill Hickey",
   "project": "DNA extracts from soil samples",
   "workdir": "/mnt/grl/brc/data/hickey_bill/rt319_soil_community/qiime_analysis_sumo/",
   "OTU_picking": "open",
   "reference": "16S",
   "data_dir": "/mnt/grl/brc/data/hickey_bill/rt319_soil_community/data",
   "delivery": "/mnt/grl/agar/users/Hickey/projects/rt319_soil_community",
   "group": {
         "1": ["1a","1b","1c"],
         "2": ["2a","2b","2c"]
   },
   "samples": {
         "1a": ["1a_GAATTCGT-CCTATCCT_L001_R1_001.fastq.gz", "1a_GAATTCGT-CCTATCCT_L001_R2_001.fastq.gz"],
         "1b": ["1b_GAATTCGT-GGCTCTGA_L001_R1_001.fastq.gz", "1b_GAATTCGT-GGCTCTGA_L001_R2_001.fastq.gz"],
         "1c": ["1c_GAATTCGT-AGGCGAAG_L001_R1_001.fastq.gz", "1c_GAATTCGT-AGGCGAAG_L001_R2_001.fastq.gz"],
         "2a": ["2a_GAATTCGT-TAATCTTA_L001_R1_001.fastq.gz", "2a_GAATTCGT-TAATCTTA_L001_R2_001.fastq.gz"],
         "2b": ["2b_GAATTCGT-CAGGACGT_L001_R1_001.fastq.gz", "2b_GAATTCGT-CAGGACGT_L001_R2_001.fastq.gz"],
         "2c": ["2c_GAATTCGT-GTACTGAC_L001_R1_001.fastq.gz", "2c_GAATTCGT-GTACTGAC_L001_R2_001.fastq.gz"]
   }, 
   "sequence_type": "paired",
   "sample_sequence_type": {"MS-A1": "paired", "MS-A2": "single"},
   "paired_min_overlap": "10",
   "paired_max_overlap": "250",
   "min_amplicon_length": "100",
   "max_amplicon_length": "480",
   "F_primer": "CTTGGTCATTTAGAGGAAGTAA",
   "R_primer": "GCTGCGTTCTTTATCGATGC",
   "read_cap": "0"
 }


The run configuaration needs to be in JSON formate. Specially pay attension to the comma if you have not used it often.

PIname, project: are simple descriptions of the project itself. the content will be copied to the 
             report.html.
workdir:     The directory for output. will be created if not exists
data_dir:    The directory of raw sequence data. I usually create a directory and create symbolic
             links to each sequence files.
reference:   need to be set for different type of project. specify 16S for bacteria, ITS for fungal, 
             or 18S for eukaryotic samples. The current pipeline will use greengene 2013-08 for 16S; 
             Unite 2014 for ITS; and Silva-111 for 16S and 18S reference database.
             use "denovo" if no reference database is available.
OTU_picking: Currently it has no effect. Default using "open" OTU_picking; unless "reference" is set to "denovo"
group:       are useful for visual differentiation in the plot. it should be following format:
             "groupName": ["sample1", "sample2"]
samples:     define the sample name and corresponding files. The filenames will be appended to the
             "data_dir"
sequence_type: can be paired, or single. This overrides sample_sequence_type 
sample_sequence_type: if it is different for each sample, list individually. (Not sure how useful this  
             option is. I have never used it since implementation)
read_cap:    is the number of reads to be used in analysis. can use it to use the same number of reads 
             for all samples. Use '0' for no cap, ie. use all the reads from the sample


Project files
=============
.
├── qiime_params.py                 <- default qiime parameters. Is imported to qiime.sm file
├── qiime_pipeline_manual.txt       <- old qiime pipeline manual. obselete
├── qiime.sm                        <- Main pipeline
├── report.html                     <- report html file template for report module
├── report.py                       <- custom report module
├── run_config.json                 <- example run configuration file
├── run_qiime_sge.sh                <- run qiime on SGE (wont work for ploting rules)
├── run_qiime.sh                    <- run qiime on local machine
├── toolsinfo                       <- default paths to programs used
└── Workflow.png                    <- workflow flowchart
