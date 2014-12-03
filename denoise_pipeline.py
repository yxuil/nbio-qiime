import os
import getopt
import sys
import time
import run_qiime
import glob
import math


script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline_jiang/pipeline_454denoising/'
sys.path.append(script_path)

class DenoisePipeline(object):
  input = None
  output = None
  map = None
  min_seq_length = None
  max_seq_length = None
  min_qual_score = 20
  qual_score_window = 50
  discard_bad_windows = None
  num_cpus = 2
  help = False
  
  def get_params(self):
    letters = 'i:,o:,m:,l:,L:,s:,w:,g:,n:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'min_seq_length=','max_seq_length=', 'min_qual_score=', 'qual_score_window=', 'discard_bad_windows=', 'num_cpus=', 'help']
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
      elif o in ['-l', '--min_seq_length']:
        if p.isdigit():
          self.min_seq_length = p
        else:
          raise Exception("Minimum length is not a number")
      elif o in ['-L', '--max_seq_length']:
        if p.isdigit():
          self.max_seq_length = p
        else:
          raise Exception("Maximum length is not a number")
      elif o in ['-s', '--min_qual_score']:
        if p.isdigit():
          self.min_qual_score = p
        else:
          raise Exception("Minimum quality score is not a number")
      elif o in ['-w', '--qual_score_window']:
        if p.isdigit():
          self.qual_score_window = p
        else:
          raise Exception("The quality score window is not a number")
      elif o in ['-g', '--discard_bad_windows']:
        self.discard_bad_windows = True
      elif o in ['-n', '--num_cpus']:
        if p.isdigit():
          self.num_cpus = p
        else:
          raise Exception("Number of CPUs is not a number")
      elif o in ['-h', '--help']:
        self.print_help()
        self.help = True
      else:
        self.print_help()
        self.help = True
        print o + " is not a valid parameter."
        raise Exception("Improper Input to converter")
    if self.help is False and (self.input is None or self.output is None or self.map is None or self.min_seq_length is None or self.max_seq_length is None):
      raise Exception("A required parameter is missing.")

     
  def run_pipeline(self):
    try:
      os.system("source /mnt/software/qiime/qiime_install/activate.sh")
      print "Running sffinfo tools..."
      os.system("mkdir " + self.output + "/sffinfo_output") 
      os.system("/mnt/grl/brc/application/qiime_pipeline_jiang/pipeline_454denoising/sffinfo " + self.input + " > " + self.output + "/sffinfo_output/sffinfo.sff.txt")
      os.system("/mnt/software/454/apps/amplicons/bin/sffinfo -s " + self.input + " > " + self.output + "/sffinfo_output/sffinfo.fasta")
      os.system("/mnt/software/454/apps/amplicons/bin/sffinfo -q " + self.input + " > " + self.output + "/sffinfo_output/sffinfo.qual")
      print "Checking map file..." 
      os.system("check_id_map.py -b -m " + self.map + " -o " + self.output +"/map_output/")
      print "Running split_libraries..." 
      os.system("split_libraries.py -b 0 -z truncate_only -f - w 50 -g" + self.output + "/sffinfo_output/sffinfo.fasta" + " -q " + self.output + "/sffinfo_output/sffinfo.qual" + " -m " + self.map + " -o " + self.output + "/split_lib_output/" + " -l " + str(self.min_seq_length) + " -L " + str(self.max_seq_length) + " -s " + str(self.min_qual_score))      
      print "Running denoise_wrapper..." 
      os.system("denoise_wrapper.py -v -i " + self.output + "/sffinfo_output/sffinfo.sff.txt -f " + self.output + "/split_lib_output/seqs.fna -m " + self.map + " -o " + self.output + "/denoise_wrapper_output/" + " -n " + str(self.num_cpus))
      print "Running inflate_denoiser..." 
      os.system("inflate_denoiser_output.py -c " + self.output + "/denoise_wrapper_output/centroids.fasta -s " + self.output + "/denoise_wrapper_output/singletons.fasta -f " + self.output + "/split_lib_output/seqs.fna -d " + self.output + "/denoise_wrapper_output/denoiser_mapping.txt -o " + self.output + "/inflate_denoiser_output/denoised.fna")
      print "Data denoising is finished."
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
  DenoisePipeline()
