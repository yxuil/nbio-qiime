import os
import getopt
import sys
import time
import glob
import math


script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline_jiang/pipeline_454denoising/'
sys.path.append(script_path)

class DenoisePipeline(object):
  input = None
  output = None
  map = None
  num_cpus = 2
  help = False
  
  def main(self):
    try:
      os.system("source /mnt/software/qiime/qiime_install/activate.sh")
      self.get_params()
      self.run_commands()
    except Exception as e:
      print e
      self.print_help()
      sys.exit()

  def get_params(self):
    letters = 'i:,o:,m:,n:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'num_cpus=', 'help']
    options, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)
    for o,p in options:
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
          raise Exception("Number of CPUs is not a number")
      elif o in ['-h', '--help']:
        self.print_help()
        self.help = True
        sys.exit()
      else:
        self.print_help()
        self.help = True
        print o + " is not a valid parameter."
        raise Exception("Improper Input to converter")
    if self.help is False and (self.input is None or self.output is None or self.map is None):
      raise Exception("A required parameter is missing.")

     
  def run_commands(self):

      sffname = self.input.split("/")[-1]
      print sffname.split(".")[0]
      print "Processing .sff files..."
      os.system("process_sff.py -f -i " + self.input + " -o " + self.output + "/process_sff_output")
       

      print "Checking map file..." 
      os.system("check_id_map.py -b -m " + self.map + " -o " + self.output +"/map_output/")

      print "Running split_libraries..." 
      os.system("split_libraries.py -o " + self.output + "/split_libraries_output -f " + self.output + "/process_sff_output/" + sffname.split(".")[0] + ".fna -q " + self.output + "/process_sff_output/" + sffname.split(".")[0] + ".qual -m " + self.map + " -b 0 -l 200 -L 400 -w 50 -g -z truncate_only")      

      print "Running denoise_wrapper..." 
      os.system("denoise_wrapper.py -v -i " + self.output + "/process_sff_output/V1_V2pool1.txt -f " + self.output + "/split_libraries_output/seqs.fna -o " + self.output + "/denoise_wrapper_output -m " + self.map)

      print "Running inflate_denoiser..." 
      os.system("inflate_denoiser_output.py -c " + self.output + "/denoise_wrapper_output/centroids.fasta -s " + self.output + "/denoise_wrapper_output/singletons.fasta -f " + self.output + "/split_libraries_output/seqs.fna -d " + self.output + "/denoise_wrapper_output/denoiser_mapping.txt -o " + self.output + "/" + sffname.split(".")[0] + "_denoisted.fna")

      print "Data denoising is finished."

  def __init__(self):
    self.main()



  def print_help(self):
    print '''
    denoise_pipeline.py runs the qiime pipeline for 16S rRNA metagenomic analysis to denoise demultiplexed 454 data. It takes a single .sff flie and outputs a denoised .fna file.
    Parameters:
    (-i, --input-dir) Required. The .sff files to be denoised.
    (-o, --output-dir) Required. The directory where you want the output data to be located.
    (-m, --map) Required. The tab-delimited qiime mapping file that contains metadata for the sample. Please read the manual for specifications on this file.
    (-n, --num_cpus) Optional. The number of cores to use for the parallel portions of the pipeline. Default is 2.
    (-h, --help) Display this help dialogue and exit.  
    '''

if(__name__ == "__main__"):
  DenoisePipeline()
