#!/mnt/software/epd/bin/ipython

import os
import getopt
import sys
import time
import glob
import csv
import collections
import math

script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline.0.1/'
sys.path.append(script_path)

class QiimeReport(object):
  input = None
  output = None
  map = None
  num_cpus = 2
  reference = None
  help = False

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
      self.run_qiime()
      self.generate_report()

  def __init__(self):
    self.main()


  def get_params(self):
    letters = 'i:,o:,m:,n:,r:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'num_cpus=', 'reference=', 'help']
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
      elif o in ['-n', '--num_cpus']:
	if p.isdigit():
	  self.num_cpus = p
	else:
	  raise Exception("num_cpus is not a number")
      elif o in ['-r', '--reference']:
	if not os.path.exists(os.path.abspath(p)):
	  raise Exception("Reference file does not exist.")
	self.reference = os.path.abspath(p)
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



  def run_qiime(self):
        
    if(self.reference is None):
      print ('de novo OTU picking...')
      os.system("pick_de_novo_otus.py -i " + self.input + " -o " + self.output + "/otus -f -a -O " + str(self.num_cpus))
      os.system("sort_otu_table.py -i " + self.output + "/otus/otu_table.biom -s SampleID -m " + self.map + " -o " + self.output +"/otus/otu_table.biom")
    else:
      print ('Open-reference OTU picking...') 
      ######the line below is added by Fenglou for debugging because it is missing -p, which seems important for fungi its analysis, use next statement for 16sRNA analysis#############
      print("pick_open_reference_otus.py -i " + self.input + " -o " + self.output + "/otus -r " + self.reference + " -p /mnt/grl/brc/application/qiime_pipeline.0.1/its_pick_open_reference_otus_params.txt -f -a -O " + str(self.num_cpus))
      os.system("pick_open_reference_otus.py -i " + self.input + " -o " + self.output + "/otus -r " + self.reference + " -p /mnt/grl/brc/application/qiime_pipeline.0.1/its_pick_open_reference_otus_params.txt -f -a -O " + str(self.num_cpus))
      ######the line above is added by Fenglou for debugging#############
      #print("pick_open_reference_otus.py -i " + self.input + " -o " + self.output + "/otus -r " + self.reference + " -f -a -O " + str(self.num_cpus))
      #os.system("pick_open_reference_otus.py -i " + self.input + " -o " + self.output + "/otus -r " + self.reference + " -f -a -O " + str(self.num_cpus))
      print("sort_otu_table.py -i " + self.output + "/otus/otu_table_mc2_w_tax.biom -s SampleID -m " + self.map + " -o " + self.output +"/otus/otu_table.biom")
      os.system("sort_otu_table.py -i " + self.output + "/otus/otu_table_mc2_w_tax.biom -s SampleID -m " + self.map + " -o " + self.output +"/otus/otu_table.biom")

    print "Creating biom_table_summary..."
    print("print_biom_table_summary.py -i " + self.output + "/otus/otu_table.biom > " + self.output + "/otus/biom_summary.txt")
    os.system("print_biom_table_summary.py -i " + self.output + "/otus/otu_table.biom > " + self.output + "/otus/biom_summary.txt")


    print "Making OTU heatmap..."
    print("make_otu_heatmap_html.py -i " + self.output + "/otus/otu_table.biom -o " + self.output + "/heatmap -n 1")
    os.system("make_otu_heatmap_html.py -i " + self.output + "/otus/otu_table.biom -o " + self.output + "/heatmap -n 1")

    print "Creating OTU network..."
    print("make_otu_network.py -m " + self.map + " -i " + self.output + "/otus/otu_table.biom -o " + self.output)
    os.system("make_otu_network.py -m " + self.map + " -i " + self.output + "/otus/otu_table.biom -o " + self.output)

    print "Creating taxa summary plots..."
    print("summarize_taxa_through_plots.py -f -i " + self.output + "/otus/otu_table.biom -o " + self.output + "/taxa_summary -m " + self.map + " -p /mnt/grl/brc/application/qiime_pipeline.0.1/summarize_taxa_params.txt")
    os.system("summarize_taxa_through_plots.py -f -i " + self.output + "/otus/otu_table.biom -o " + self.output + "/taxa_summary -m " + self.map + " -p /mnt/grl/brc/application/qiime_pipeline.0.1/summarize_taxa_params.txt")
 
    print "Creating rarefaction plots..."
    print("alpha_rarefaction.py -f -i " + self.output + "/otus/otu_table.biom -m " + self.map + " -o " + self.output +
      "/arare/ -p /mnt/grl/brc/application/qiime_pipeline.0.1/alpha_params.txt -t " + self.output + "/otus/rep_set.tre -a -O " +str(self.num_cpus))
    os.system("alpha_rarefaction.py -f -i " + self.output + "/otus/otu_table.biom -m " + self.map + " -o " + self.output +
      "/arare/ -p /mnt/grl/brc/application/qiime_pipeline.0.1/alpha_params.txt -t " + self.output + "/otus/rep_set.tre -a -O " +str(self.num_cpus))

    print "Creating beta diversity plots..."
    # use the minimum reads number as the parameter e for even sampling
    for item in open (self.output + "/otus/biom_summary.txt"):
      if "Min" in item:
        self.min_depth =str(int(float(item.strip()[5:])))
    print ('Depth of coverage for even sampling is ' + self.min_depth)   
    print("beta_diversity_through_plots.py -f -i " + self.output + "/otus/otu_table.biom -m " + self.map + " -o " + self.output + "/bdiv -t " + self.output + "/otus/rep_set.tre -e " + self.min_depth)
    os.system("beta_diversity_through_plots.py -f -i " + self.output + "/otus/otu_table.biom -m " + self.map + " -o " + self.output + "/bdiv -t " + self.output + "/otus/rep_set.tre -e " + self.min_depth)


  def print_help(self):
    print '''
    qiime_report.py is a part of the qiime pipeline for 16S rRNA metagenomic analysis. It can also run by itself using a single FNA file as input.
    Required Parameters:
    (-i, --input-dir) The input FNA file.
    (-o, --output-dir) The output directory where you want the results to be located.
    (-m, --map) The tab-delimited qiime mapping file that contains metadata for the samples to be analyzed.
    Optional Parameters:
    (-r, --reference) Reference database for OTU clustering. It should be in fasta format. If a reference is given, pick_open_reference_otus.py runs. If not, pick_de_novo_otus.py runs.
    (-n, --num_cpus) The number of cores to use for the parallel portions of the pipeline. [default: 2].
    Other Options: 
    (-h, --help) Display this help dialogue and exit.
    For complete information about this pipeline, please refer to the manual.
    '''
 
  def generate_report(self):
    outputDir= os.path.abspath(self.output)
    PIName = "PI Name"
    ref_db = "greengenes"

    print "Generating report..."
    os.system("cp /mnt/grl/brc/application/qiime_pipeline.0.1/Workflow.png " + outputDir + "/Workflow.png")
    f=open(outputDir+"/report.html","w")
    f.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n")
    f.write("<html>\n")
    f.write("<head>\n")
    f.write("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n")

    f.write("<style type=\"text/css\">\n")
    f.write("body, th, td, p { font-family: helvetica, sans-serif; font-size: small }\n")
    f.write("h1 { font-size: medium; text-align: center }\n")
    f.write("h2, h3 { font-size: small; text-align: left }\n")
    f.write("table {}\n") 
     
    f.write(".brc { border: 1px solid black; border-collapse: collapse; margin: 10px 10px 15px 0; width: 820px }")


    f.write("th { text-align: left; color: white; background-color: #aa3333; padding: 3px 15px 3px 3px; border-bottom: 1px solid black; width: 150px }\n")
    f.write("td { padding: 3px 3px 3px 4px; border-bottom: 1px solid gray  }\n")
    f.write(".ntitle { color: black; font-family:Arial,Verdana; font-size:11; font-weight:bold;}")
    f.write(".left { text-align: left }\n")
    f.write(".right { text-align: right }\n")
    f.write(".center { text-align: center }\n")
    f.write(".boldcentered { text-align: center; font-weight: bold }\n")
    f.write("</style>\n")
    

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

# scan through the map file and create 'sample' objects that contain all the metadata about each sample. 
def parse_map(map):

  # Find the indexes of all metadata that doesn't have a definite location. 
  fprimer_index = 2 
  rprimer_index = -1
  treatment_index = -1
  filename_index = -1
  sampleID_index = 0 
  barcode_index = 1 
  description_index = -1
  other_indexes = []
  sample_list =[] 
  Sample = collections.namedtuple('Sample', ['fprimer', 'rprimer', 'treatment', 'filename', 'sampleID', 'barcode', 'description', 'other'])
  with open(map, "rb") as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    first_line = reader.next()
    # read header line and figure out the indexes of the metadata
    for i in range(len(first_line)):
      if first_line[i].lower() == "treatment":
        treatment_index = i 
      elif first_line[i].lower() == "reverseprimer":
        rprimer_index = i 
      elif first_line[i].lower() == "filename":
        filename_index = i 
      elif first_line[i].lower() == "description":
        description_index = i 
      else:
        if i!=fprimer_index and i!=sampleID_index and i!=barcode_index:
          other_indexes += [i] 

    # go through the rows of the map and generate Sample objects for each sample
    for row in reader:
      if row[0][0] != '#':
        sample_data = Sample(fprimer = row[fprimer_index], rprimer = row[rprimer_index], treatment = row[treatment_index], filename = row[filename_index],
          sampleID = row[sampleID_index], barcode = row[barcode_index], description = row[description_index], other = []) 
        for i in other_indexes:
          sample_data = sample_data._replace(other= sample_data.other + [row[i]])
        sample_list += [sample_data]
  # return a list containing the metadata for all the samples. 
  return sample_list



if(__name__ == "__main__"):
  QiimeReport()

