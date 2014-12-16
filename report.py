def generate_report(config):
    outputDir= config["workdir"]
    PIName = "PI Name"
    ref_db = "greengenes"

    print "Generating report..."
    os.system("cp /mnt/grl/brc/application/qiime_pipeline.0.1/Workflow.png " + outputDir + "/Workflow.png")
    f=open(outputDir+"/report.html","w")
    f.write("""<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">
    <html>
    <head>
    <meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">

    <style type=\"text/css\">
    body, th, td, p { font-family: helvetica, sans-serif; font-size: small }
    h1 { font-size: medium; text-align: center }
    h2, h3 { font-size: small; text-align: left }
    table {}

    .brc { border: 1px solid black; border-collapse: collapse; margin: 10px 10px 15px 0; width: 820px }")


    th { text-align: left; color: white; background-color: #aa3333; padding: 3px 15px 3px 3px; border-bottom: 1px solid black; width: 150px }
    td { padding: 3px 3px 3px 4px; border-bottom: 1px solid gray  }
    .ntitle { color: black; font-family:Arial,Verdana; font-size:11; font-weight:bold;}")
    .left { text-align: left }
    .right { text-align: right }
    .center { text-align: center }
    .boldcentered { text-align: center; font-weight: bold }
    </style>


    </head>
    <body>
    <div class=\"boldcentered\">UW Biotechnology Center Bioinformatics Resource Center - "+time.asctime()+"</div>
    f.write('<div style="margin-top: 10px" class="center">Please direct questions to: <strong>brc@biotech.wisc.edu</strong></div>'+"
    <hr>
    <ol>
    <li><a href='#summary'>Summary</a></li>
    <li><a href='#workflow'>Workflow</a></li>
    <li><a href='#taxa_summary'>Taxa Summary</a></li>
    <li><a href='#arare'>Alpha Rarefaction</a></li>
    <li><a href='#bdiv'>Beta Diversity</a></li>
    <li><a href='#phylogeny'>Phylogenetic Tree</a></li>
    <li><a href='#heatmap'>OTU Heatmap</a></li>
    <li><a href='#software'>Software</a></li>
    </ol>
    <hr>
    """)



    samples = parse_map(self.map)
    #summary
    f.write("<div class=\"left\"><strong><a id=summary name=summary>Summary</a></strong></div>")
    f.write("<table class=\"brc\">\n")
    #f.write("<caption>Summary</caption>\n")
    f.write("<tr><th scope='row'>PI Name</th><td>"+PIName+"</td></tr>\n")
    f.write("<tr><th scope='row'>Project Description</th><td>"+samples[0].description+"</td></tr>\n")
    sample_names = ""
    for s in samples:
      sample_names += s.sampleID + ", "
    f.write("<tr><th scope='row'>Samples</th><td>"+sample_names+"</td></tr>\n")
    if(not (self.reference is None)):
      f.write("<tr><th scope='row'>Reference Database</th><td>"+ self.reference +"</td></tr>\n")
    #f.write("<tr><th>Report Generation Date</th><td>"+time.asctime()+"</td></tr>\n")
    f.write("</table>\n")
    f.write("<hr>\n");
    #reads and coverage
    '''
    f.write("<div class='left'><strong><a id='reads' name='reads'>Reads statistics</a></strong> (click the links in mapped reads field to download the final bam files)</div>")
    f.write("<table>\n")
    f.write("<tr>\n")
    f.write("<th scope='col'>Sample Name</th>\n")
    f.write("<th scope='col'>Total Reads</th>\n")
    f.write("<th scope='col'>Mapped Reads</th>\n")
    f.write("<th scope='col'>Mapped Reads Percentage</th>\n")
    f.write("<th scope='col'>Coverage</th>\n")
    '''
    # So I need to figure out where to find the read statistics.


    # The workflow
    f.write("<div class='left'><strong><a id='workflow' name='workflow'>Workflow Chart</a></strong></div><table>\n")
    f.write('<tr>\n')
    f.write('<td class=center><img src="Workflow.png" alt="Workflow Chart"></td></tr>\n')
    f.write('</table>\n')
    f.write('<hr>\n')


    path = ('./taxa_summary/taxa_summary_plots/bar_charts.html')
    f.write("<div class=\"left\"><strong><a id=taxa_summary name=taxa_summary>Taxa Summary</a></strong></div>")
    f.write("<li><a href=\'" + path + "\'>View Full Interactive Taxa Summary</a></li>\n");

    with open(outputDir + '/taxa_summary/taxa_summary_plots/bar_charts.html', 'r') as bar_summary:
      start=False
      for line in bar_summary:
        if start:
  	  str = line.replace('a href=\'', 'a href=\'' + 'taxa_summary/taxa_summary_plots/')
	  str = str.replace('a href=\"', 'a href=\"' + 'taxa_summary/taxa_summary_plots/')
	  str = str.replace('img src=\'', 'img src=\'' + 'taxa_summary/taxa_summary_plots/')
	  f.write(str)
        if line.find("otu_table_L6.txt") != -1:
	  start = True

    f.write("<hr>\n");


    f.write("<div class=\"left\"><strong><a id=arare name=arare>Alpha Rarefaction</a></strong></div>")
    path = ('./arare/alpha_rarefaction_plots/rarefaction_plots.html')
    f.write("<li><a href=\'" + path + "\'>View Alpha Rarefaction Plots</a></li>\n");

    f.write("<hr>\n")
    f.write("<div class=\"left\"><strong><a id=bdiv name=bdiv>Beta Diversity</a></strong></div>")


    path = ('./bdiv/unweighted_unifrac_2d_continuous/unweighted_unifrac_pc_2D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 2D Continuous Unweighted Unifrac </a></li>\n");

    path = ('./bdiv/unweighted_unifrac_2d_discrete/unweighted_unifrac_pc_2D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 2D Discrete Unweighted Unifrac </a></li>\n");

    path = ('./bdiv/weighted_unifrac_2d_continuous/weighted_unifrac_pc_2D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 2D Continuous Weighted Unifrac </a></li>\n");

    path = ('./bdiv/weighted_unifrac_2d_discrete/weighted_unifrac_pc_2D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 2D Discrete Weighted Unifrac </a></li>\n");

    path = ('./bdiv/unweighted_unifrac_3d_continuous/unweighted_unifrac_pc_3D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 3D Continuous Unweighted Unifrac</a></li>\n");

    path =  ('./bdiv/unweighted_unifrac_3d_discrete/unweighted_unifrac_pc_3D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 3D Discrete Unweighted Unifrac </a></li>\n");

    path = ('./bdiv/weighted_unifrac_3d_continuous/weighted_unifrac_pc_3D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 3D Continuous Weighted Unifrac </a></li>\n");

    path = ('./bdiv/weighted_unifrac_3d_discrete/weighted_unifrac_pc_3D_PCoA_plots.html')
    f.write("<li><a href=\'" + path + "\'>View 3D Discrete Weighted Unifrac </a></li>\n");

    f.write("<hr>\n");


    #Heatmap
    f.write("<div class='left'><strong><a id='heatmap' name='heatmap'>Heatmap</a></strong></div>")
    path = './heatmap/otu_table.html'
    f.write("<li><a href=\'" + path + "\'>View OTU Heatmap</a></li>\n");
    f.write("<hr>\n")


    # Phylogeny
    f.write("<div class='left'><strong><a id='phylogeny' name='phylogeny'>Phylogenetic Tree</a></strong></div>")
    path = ('./otus/rep_set.tre')
    f.write("<li><a href=\'" + path + "\'>Right click and select \"save as\" to download tree</a></li>\n");
    f.write("<hr>\n")



    #software
    f.write("<div class='left'><strong><a id='software' name='software'>Software</a></strong></div>")
    f.write("<table class=\"brc\">\n")
    f.write("<tr><th>Main Pipeline</th><td>QIIME</td></tr>\n")
    f.write("<tr><th>OTU Clustering</th><td>Uclust</td></tr>\n")
    f.write("<tr><th>Quality Filtering</th><td>Usearch</td></tr>\n")
    f.write("<tr><th>Alignment</th><td>PyNAST</td></tr>\n")
    f.write("<tr><th>Taxonmy Assignment</th><td>RDP Classifier</td></tr>\n")
    f.write("<tr><th>Phylogeny Generation</th><td>FastTree</td></tr>\n")
    f.write("</table>\n")
    f.write("<hr>\n");


    print "Report Finished"


if __name__ == "__main__":
    import sys, json
    config = json.load(sys.argv[1])
    generate_report(config)