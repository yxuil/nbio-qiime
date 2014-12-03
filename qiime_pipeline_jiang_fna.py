
import os
import getopt
import sys
import time
import run_qiime
import glob

script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline_jiang/pipeline_454denoising/'
sys.path.append(script_path)

class Pipeline(object):
  input = None
  output = None
  map = None
  cores = 2
  reference = None
  help = False
  sort = 'SampleID'

  def get_params(self):
    letters = 'i:,o:,m:,c:,r:,s:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'cores=', 'reference=', 'sort=', 'help']
    opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)
    for o,p in opts:
      if o in ['-i', '--input-dir']:
        if not os.path.exists(os.path.abspath(p)):
	  raise Exception("Invalid input path")
	self.input = os.path.abspath(p)
      elif o in ['-o', '--output-dir']:
	if not os.path.exists(os.path.abspath(p)):
	  os.makedirs(os.path.abspath(p))
	self.output = os.path.abspath(p)
      elif o in ['-m', '--map']:
	if not os.path.exists(os.path.abspath(p)):
	  raise Exception("Map file does not exit")
	self.map = os.path.abspath(p)
      elif o in ['-c', '--cores']:
	if p.isdigit():
	  self.cores = p
	else:
	  raise Exception("Cores is not a number")
      elif o in ['-r', '--reference']:
	if not os.path.exists(os.path.abspath(p)):
	  raise Exception("Reference file does not exist.")
	self.reference = os.path.abspath(p)
      elif o in ['-s', '--sort']:
	self.sort = p
      elif o in ['-h', '--help']:
	self.print_help()
	self.help = True
      else:
	self.print_help()
	self.help = True
	print o + " is not a valid parameter."
	raise Exception("Improper Input to converter")
    if self.help is False and (self.input is None or self.output is None or self.map is None):
      raise Exception("A required parameter is missing!")

    # if I run these by creating objects of them, I can handle errors that they produce. 
  def run_pipeline(self):
    try:
      run_qiime_string = os.path.join(script_path, "run_qiime.py") + " -- -i " + self.input + " -o " + self.output + " -m " + self.map + " -c " + str(self.cores) + " -s " + self.sort
      if (not (self.reference is None)):
        run_qiime_string += " -r " + self.reference
      os.system(run_qiime_string)
      self.generate_HTML()
    except Exception as e:
      print e
      sys.exit(1)

  def main(self):
    try:
      self.get_params()
    except Exception as e:
      print e
      self.print_help()
      print "\n\n"
      print e
      sys.exit()
    if(not self.help):
      self.run_pipeline()

  def __init__(self):
    self.main()


    # generate an HTML summary report for the qiime analysis. 
    #def generate_report():






  def print_help(self):
    print '''
    qiime_pipeline.py runs the qiime pipeline for 16S rRNA metagenomic analysis on given input data. 
    Parameters:
    (-i, --input-dir) The input .fna file to be analyzed. The .fna should be demultiplexed, denoised and filtered by quality score.
    (-o, --output-dir) The directory where you want the output data and report to be located.
    (-m, --map) The tab-delimited qiime mapping file that contains metadata for the samples. Please read the manual for specifications on this file
      The manual is located in ???
    (-c, --cores) Optional. The number of cores to use for the parallel portions of this pipeline. Default is 2. 
    (-r, --reference) Optional. A reference database (such as greengenes) for OTU clustering. If none is specified, de novo clustering will be used.
    (-s, --sort) Optional. The metadata parameter by which to sort the samples for output. Defaults to sampleID. 
    (-h, --help) Display this help dialogue and exit. 
    
    If this is your first time using the pipeline, please read the manual. It will save you an incredible amount of headache
    trying to format the map file. 
    '''

  # Generate the HTML report from the pipeline. 
  def generate_HTML(self):
    outputDir= os.path.abspath(self.output)
    PIName = "PI Name"
    ref_db = "greengenes"


    os.system("mkdir -p "+outputDir+"/report_files/")
    os.system("cp /mnt/software/qiime/qiime_pipeline_files/Workflow.png " + outputDir + "/report_files/Workflow.png")
    f=open(outputDir+"/report.html","w")
    f.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n")
    f.write("<html>\n")
    f.write("<head>\n")
    f.write("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n")

    f.write("<style type=\"text/css\">\n")
    f.write("body, th, td, p { font-family: helvetica, sans-serif; font-size: small }\n")
    f.write("h1 { font-size: medium; text-align: center }\n")
    f.write("h2, h3 { font-size: small; text-align: left }\n")
    f.write("table {}\n") # I could change this to get rid of the table border. make two types of table, one with a border and the default
    # without. then add tags to all the tables I use so that they have the borders and then the qiime tables will follow the default table behavior. 
    f.write(".brc { border: 1px solid black; border-collapse: collapse; margin: 10px 10px 15px 0; width: 820px }")

    #f.write("table tr:nth-child(even) td { background-color: #eeeeee }\n")# uncommenting this fixes it.
    f.write("th { text-align: left; color: white; background-color: #aa3333; padding: 3px 15px 3px 3px; border-bottom: 1px solid black; width: 150px }\n")
    f.write("td { padding: 3px 3px 3px 4px; border-bottom: 1px solid gray  }\n")
    f.write(".ntitle { color: black; font-family:Arial,Verdana; font-size:11; font-weight:bold;}")
    f.write(".left { text-align: left }\n")
    f.write(".right { text-align: right }\n")
    f.write(".center { text-align: center }\n")
    f.write(".boldcentered { text-align: center; font-weight: bold }\n")
    f.write("</style>\n")
    #f.write("<link rel=\"stylesheet\" href=\"./css/qiime_style.css\" type=\"text/css\">")




    f.write("</head>\n")
    f.write("<body>\n")
    f.write("<div class=\"boldcentered\">UW Biotechnology Center Bioinformatics Resource Center - "+time.asctime()+"</div>\n")
    f.write('<div style="margin-top: 10px" class="center">Please direct questions to: <strong>brc@biotech.wisc.edu</strong></div>'+"\n");
    f.write("<hr>\n");
    f.write("<ol>\n");
    f.write("<li><a href='#summary'>Summary</a></li>\n");
    f.write("<li><a href='#workflow'>Workflow</a></li>\n");
    f.write("<li><a href='#taxa_summary'>Taxa Summary</a></li>\n");
    f.write("<li><a href='#arare'>Alpha Rarefaction</a></li>\n");
    f.write("<li><a href='#bdiv'>Beta Diversity</a></li>\n");
    f.write("<li><a href='#phylogeny'>Phylogenetic Tree</a></li>\n")
    f.write("<li><a href='#heatmap'>OTU Heatmap</a></li>\n")
    f.write("<li><a href='#software'>Software</a></li>\n");
    f.write("</ol>\n");
    f.write("<hr>\n");



    samples = run_qiime.parse_map(self.map)
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
    f.write('<td class=center><img src="report_files/Workflow.png" alt="Workflow Chart"></td></tr>\n')
    f.write('</table>\n')
    f.write('<hr>\n')
    # The taxonomy summary
    # How about I put a table of all the stuff in the excell sheet and include the pictures. Basically just replicate the html page that the pipeline creates. 
    # I could just link to it but it would look much worse, let's see how much work it'd be to copy over the stuff from the html that qiime generates. 
    # okay. I think I can just hard code the header and scripts from the qiime output then just grab everything after the line containing "otu_table_L5.txt"

    # how do I make this happen? 
    # I'll have to use a filereader and iterate though line by line until I find the string I'm looking for
    # or there might be a find method that returns line numbers. 

    path = ('./otus/taxa_summary/taxa_summary_plots/bar_charts.html')
    f.write("<div class=\"left\"><strong><a id=taxa_summary name=taxa_summary>Taxa Summary</a></strong></div>")
    f.write("<li><a href=\'" + path + "\'>View Full Interactive Taxa Summary</a></li>\n");

    with open(outputDir + '/otus/taxa_summary/taxa_summary_plots/bar_charts.html', 'r') as bar_summary:
      start=False
      for line in bar_summary:
        if start:
  	  str = line.replace('a href=\'', 'a href=\'' + 'otus/taxa_summary/taxa_summary_plots/')
	  str = str.replace('a href=\"', 'a href=\"' + 'otus/taxa_summary/taxa_summary_plots/')
	  str = str.replace('img src=\'', 'img src=\'' + 'otus/taxa_summary/taxa_summary_plots/')
	  f.write(str)
        if line.find("otu_table_L5.txt") != -1:
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
    path = './otus/heatmap/otu_table.html'
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




if(__name__ == "__main__"):
  Pipeline()

