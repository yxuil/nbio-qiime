
import os
import getopt
import sys
import time
import run_qiime
import glob

script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline_jiang/'
sys.path.append(script_path)

class PGM_Pipeline(object):
  input = None
  output = None
  map = None
  num_cpus = 2
  length = None
  quality = None
  reference = None
  help = False
  sort = 'SampleID'

  def get_params(self):
    letters = 'i:,o:,m:,n:,l:,q:,r:,s:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'num_cpus=','length=', 'quality=', 'reference=', 'sort=', 'help']
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
	  raise Exception("The number of cpus is not a number")
      elif o in ['-q', '--quality']:
	if p.isdigit():
	  self.quality = p
	else:
	  raise Exception("Quality is not a number")
      elif o in ['-l', '--length']:
	if p.isdigit():
	  self.length = p
	else:
	  raise Exception("Length is not a number")
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
    if self.help is False and (self.input is None or self.output is None or self.map is None or self.quality is None or self.length is None):
	raise Exception("A required parameter is missing.")

     
  def run_trim(self):
    try:
      os.system(os.path.join(script_path, "bam_to_fastq.py") + " -- -i " + self.input + " -o " + self.output + " -m " + self.map)
      os.system(os.path.join(script_path, "trim_primers.py") + " -- -i " + self.output + "/fastq/ -o " + self.output + "/trimmed_fastq/ -l " + self.length + " -q " + self.quality + " -m " + self.map)
    except Exception as e:
      print e
      sys.exit(1)

  def run_qiime_report(self):
    run_QiimeReport_command = "python " + os.path.join(script_path, "qiime_report_j2.py") + " -i " + self.output + "/trimmed_fastq/all_reads.fna -o " + self.output + " -m " + self.map + " -n " + str(self.num_cpus)
    os.system(run_QiimeReport_command)

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
      self.run_trim()
      self.run_qiime_report()

  def __init__(self):
    self.main()



  def print_help(self):
    print '''
    qiime_pipeline.py runs the qiime pipeline for 16S rRNA metagenomic analysis on given input data. 
    Parameters:
    (-i, --input-dir) Required. The directory containing all the .bam or .fastq files to be analyzed.
    (-o, --output-dir) Required. The directory where you want the output data and report to be located.
    (-m, --map) Required. The tab-delimited qiime mapping file that contains metadata for the samples. Please read the manual for specifications on this file.
    (-q, --quality) Require. The quality threshold for reads. All reads below this quality will be removed before analysis is performed.
    (-l , --length) Required. the minimum length for reads. All reads shorter than this length will be removed before analysis is performed. 
    (-c, --cores) Optional. The number of cores to use for the parallel portions of the pipeline. Default is 2.
    (-r, --reference) Optional. A reference database (such as greengenes) for OTU clustering. If none is specified, de novo clustering will be used.
    (-s, --sort) Optional. The metadata parameter by which to sort the samples for output. Defaults to sampleID. 
    (-h, --help) Display this help dialogue and exit. 
    The complete manual is located in /mnt/software/qiime_andrew/my_script/qiime_pipeline/documentation.
    If this is your first time using the pipeline, please read the manual. It will save you an incredible amount of effort to format the map file. 
    '''


if(__name__ == "__main__"):
  PGM_Pipeline()

