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
  min_seq_length = 200
  max_seq_length = 1000
  min_qual_score = 20
  qual_score_window = 50
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
    letters = 'i:,o:,m:,n:,l:,L:,s:,w:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'num_cpus=', 'min_seq_length=', 'max_seq_length=', 'min_qual_score=', 'qual_score_window=', 'help']
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
      print ("Processing " + sffname +" ...")
      os.system("process_sff.py -f -i " + self.input + " -o " + self.output + "/process_sff_output")
       

      print "Checking map file..." 
      os.system("check_id_map.py -b -m " + self.map + " -o " + self.output +"/map_output/")

      print "Running split_libraries..." 
      os.system("split_libraries.py -b 0 -z truncate_only -g -o " + self.output + "/split_libraries_output -f " + self.output + "/process_sff_output/" + sffname.split(".")[0] + ".fna -q " + self.output + "/process_sff_output/" + sffname.split(".")[0] + ".qual -m " + self.map + " -l " + str(self.min_seq_length) + " -L " + str(self.max_seq_length) + " -s " + str(self.min_qual_score) + " -w " + str(self.qual_score_window))      

      print "Running denoise_wrapper..." 
      os.system("denoise_wrapper.py -v -i " + self.output + "/process_sff_output/" + sffname.split(".")[0] + ".txt -f " + self.output + "/split_libraries_output/seqs.fna -o " + self.output + "/denoise_wrapper_output -m " + self.map)

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
    (-m, --map) Required. The tab-delimited qiime mapping file that contains metadata for the sample. The BarcodeSequence column should be empty.
    (-n, --num_cpus) Optional. The number of cores to use for the parallel portions of the pipeline. [default: 2].
    (-l, --min-seq-length) Optional. Minimum sequence length, in nucleotides. [default: 200]
    (-L, --max-seq-length) Optional. Maximum sequence length, in nucleotides. [default: 1000]
    (-s, --min-qual-score) Optional. Min average qual score allowed in reads. [default: 20]
    (-w, --qual_score_window) Optional. Enable sliding window test of quality scores. If the average score of a continuous set of w nucleotides falls below the threshold (see -s for default), the sequence is discarded. A good value would be 50. 0 (zero) means no filtering. Must pass a .qual file (see -q parameter) if this functionality is enabled. The behavior for this function is to discard any sequences where a bad window is found. [default: 50] 
    (-h, --help) Display this help dialogue and exit.  
    '''

if(__name__ == "__main__"):
  DenoisePipeline()
