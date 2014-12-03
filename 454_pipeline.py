#!/mnt/software/epd/bin/ipython

import os
import getopt
import sys
import time
import glob
import math

script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline_jiang/'
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
      
      # The while loop is to denoise the .sff files specified in the map file one by one.     
      f1 = open(self.map)
      lines=f1.readlines()
      n=2                                       #Samples start on 3rd line
      while n < len(lines):
        self.sffname = lines[n].split('\t')[4]  #Filename is 5th column.
        self.sffid = self.sffname.split('.')[0]
        print ('The following sff file is being denoised: ' + self.sffname)
        print "  Getting metadata from the mapping file..."
        self.singlemap = self.sffid + '_map.txt'
        f2 = open (self.output + "/" + self.singlemap, 'w')
        f2.write(lines[0])                       #Header is on 1st line
        f2.write(lines[n])    
        f2.close()
        self.run_denoise()
        n+=1
      f1.close()
      print "Data denoising is finished."

      print "Combining denoised.fna files..."
      os.system("cat " + self.output + "/*denoisted.fna > " + self.output + "/combined_denoised.fna")
      
      print "Running qiime_report..."
      self.run_QiimeReport()
      
      
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



 
  def run_denoise(self):

      print ("  Processing the .sff file...")
      os.system("process_sff.py -f -i " + self.input + "/" + self.sffname + " -o " + self.output + "/" + self.sffid + "_process_sff_output")

      print "  Running split_libraries..." 
      split_lib_command = "split_libraries.py -b 0 -z truncate_only -g -o " + self.output + "/" + self.sffid + "_split_libraries_output -f " + self.output + "/" + self.sffid + "_process_sff_output/" + self.sffid + ".fna -q " + self.output + "/" + self.sffid + "_process_sff_output/" + self.sffid + ".qual -m " + self.output + "/" + self.singlemap + " -l " + str(self.min_seq_length) + " -L " + str(self.max_seq_length) + " -s " + str(self.min_qual_score) + " -w " + str(self.qual_score_window)
      print split_lib_command
      os.system(split_lib_command) 

      print "  Running denoise_wrapper..."
      denoise_wrapper_command = "denoise_wrapper.py -v -i " + self.output + "/" + self.sffid + "_process_sff_output/" + self.sffid + ".txt -f " + self.output + "/" + self.sffid + "_split_libraries_output/seqs.fna -o " + self.output + "/" + self.sffid + "_denoise_wrapper_output -m " + self.output + "/" + self.singlemap 
      print denoise_wrapper_command
      os.system(denoise_wrapper_command)

      print "  Running inflate_denoiser..."
      inflate_denoiser_command = "inflate_denoiser_output.py -c " + self.output + "/" + self.sffid + "_denoise_wrapper_output/centroids.fasta -s " + self.output + "/" + self.sffid + "_denoise_wrapper_output/singletons.fasta -f " + self.output + "/" + self.sffid + "_split_libraries_output/seqs.fna -d " + self.output + "/" + self.sffid + "_denoise_wrapper_output/denoiser_mapping.txt -o " + self.output + "/" + self.sffid + "_denoisted.fna"
      print inflate_denoiser_command
      os.system(inflate_denoiser_command)

  def run_QiimeReport(self):
    run_QiimeReport_command = "python " + os.path.join(script_path, "qiime_report_j1.py") + " -i " + self.output + "/combined_denoised.fna -o " + self.output + " -m " + self.map + " -n " + str(self.num_cpus)
    print run_QiimeReport_command
    os.system(run_QiimeReport_command)
         
       

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
