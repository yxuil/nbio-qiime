#!/mnt/software/epd/bin/ipython

import os
import getopt
import sys
import time
import glob
import math
import csv
import collections
from Bio import SeqIO


script_path = os.path.dirname(os.path.realpath(__file__)) #'/mnt/grl/brc/application/qiime_pipeline.0.1/'
sys.path.append(script_path)



class make_groups(object):
  input = None
  output = None
  group = None
  help = False




  def main(self):
    try:
      self.get_params()
      self.write_groups()

      
    except Exception as e:
      print e
      self.print_help()
      sys.exit()





  def get_params(self):
    letters = 'i:,o:,g:,h'
    keywords = ['input-dir=', 'output-dir=', 'group=', 'help']
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

      elif o in ['-g', '--group-name']:
          self.num_cpus = p
      elif o in ['-h', '--help']:
        self.print_help()
        self.help = True
        sys.exit()
      else:
        self.print_help()
        self.help = True
        print o + " is not a valid parameter."
        raise Exception("Improper Input to converter")
    if self.help is False and (self.input is None or self.output is None or self.group is None):
      raise Exception("A required parameter is missing.")


  def write_group(self)
    f1=open(self.input)
    f2=open(self.output, 'w')
    for line in f1:
      if line[0]=='>':
      seq_name = line[1:15]
      f2.write(seq_name + '\t' + self.group + '\n')
    f1.close()
    f2.close()


  def print_help(self):
    print '''
    qiime_pipeline.py runs the qiime pipeline for 16S rRNA metagenomic analysis. It takes demultiplexed SFF, BAM, and FASTQ files from 454, pgm and illumina sequencers.
    Required Parameters:
    (-i, --input-dir) The input directory containing all the data you want to analyzed. 
    (-o, --output-dir) The output directory where you want the results to be located.
    (-g, --group-name) The tab-delimited qiime mapping file that contains metadata for the samples to be analyzed.
    '''




if(__name__ == "__main__"):
  make_groups()
