#parameters for QIIME pipelines

pick_otus:enable_rev_strand_match True
pick_otus:max_accepts 1
pick_otus:max_rejects 8
pick_otus:stepwords 8
pick_otus:word_length 8

assign_taxonomy:id_to_taxonomy_fp {QIIME}/{REF}/taxonomy/{HOMOLOGY}_otu_taxonomy.txt
assign_taxonomy:reference_seqs_fp {QIIME}/{REF}/rep_set/{HOMOLOGY}_otus.fasta

summarize_taxa:level    2,3,4,5,6,7
plot_taxa_summary:chart_type    bar

alpha_diversity:metrics {aMETRICS}
#observed_species,shannon,PD_whole_tree,chao1 <- 16S
#observed_species,chao1,shannon  <- ITS

beta_diversity:metrics {bMETRICS}
#unweighted_unifrac,weighted_unifrac <- 16S
#bray_curtis <- ITS
