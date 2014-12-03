#!/mnt/software/epd/bin/ipython


#TODO:
#get it working if the user only supplies a relative path. 
#make a help dialogue and make it pop up if the user supplies bad input or the -h keyword


import collections
import getopt
import sys
import os
import csv
import run_qiime
from Bio import SeqIO #uncomment this before running
import shutil

class Trimmer(object):
  input= None
  output= None
  length= None
  quality= None
  map= None
  files=[]
  samples= []
  trimmed_files = []

  # main class just calls the two methods to process the data
  def main(self):
    os.system("source /mnt/software/qiime/qiime_install/activate.sh")
    try:
      self.get_params()
      self.samples = run_qiime.parse_map(self.map)
      self.trim()
      self.merge_files()
    except Exception as e:
      #print_help() # Uncomment if running the trimmer as a standalone app. keep commented when in pipeline
      print e
      shutil.rmtree(self.output, True)
      exit()

# what are my inputs:
# a directory of fastq files (-i, --input-dir)
# a destination directory for trimmed fastq files (-o, --output-dir)
# a minimum length (-l, --length)
# a quality threshold (-q, --quality)
  def get_params(self):
    letters = 'i:,o:,l:,q:,m:' # letters for passing args
    keywords = ['input-dir=', 'output-dir=', 'length=', 'quality=', 'map='] # keywords for passing args

    opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords) # start at second element because first is script name. opts, extraparams are lists of tuples
    for o,p in opts:
      if o in ['-i', '--input-dir']:
        self.input = p
      elif o in ['-o', '--output-dir']:
        self.output = p
      elif o in ['-l', '--length']:
        self.length = p
      elif o in ['-q', '--quality']:
        self.quality = p
      elif o in ['-m', '--map']:
        self.map = p
      else:
        raise Exception("Improper Input")

    #sanity check
    if(self.input==None or self.output==None or self.length==None or self.quality==None or self.map==None):
      raise Exception("Imporper Input")

  # just your basic init function
  def __init__(self):
    self.main()

  '''
  Trim the adaptors off the sequences 

  I need to check my paths. Make sure these guys work if the user only specifies a relative path
  '''
  def trim(self):

    #check that the specified path exists path.exists
    if (not os.path.exists(self.input)):
      raise Exception("Invalid Input Path")

    # check that the destination directory exists, if not, make it
    if (not os.path.exists(self.output)):
      os.mkdir(self.output)
      print "Directory " + self.output + " does not exist. Creating it now."

    # iterate through fq files calling cutadapt
    # I can parallelize this for loop to speed this up. 
    for s in self.samples:
      # create a log file
      log = self.output + s.sampleID + ".log"
      (open(log, 'w')).close() # do I even need to do this?????

      # create a file to hold the fq after the 5' primer is removed
      tmpfq = self.output + s.sampleID + "_tmp.fq"
      (open(tmpfq, 'w')).close()

      # remove the 5' primer
      os.system("cutadapt -g ^" + s.fprimer + " -q " + self.quality + " -o " + tmpfq + " " + self.input + s.sampleID + ".fq" + " > " + log)
      print "5\' primer removed from " + s.sampleID


      # remove the 3' primer
      os.system("cutadapt -a " + s.rprimer + " -m " + self.length + " -o " + self.output + s.sampleID + "_trimmed.fq " + tmpfq + " >> " + log)
      print "3\' primer removed from " + s.sampleID


      # destroy the file that only had the 5' primer removed
      os.system("rm " + tmpfq)



  # Go through our fastq files and merge them all into one .fna file with the qiime formatting
  def merge_files(self):
    print "merging files"
    fna_out=open(self.output + "/all_reads.fna", 'w') # it would be nice to generalize the name

    for s in self.samples:
      print "merging " + s.sampleID
      f_in = open(self.output + s.sampleID + "_trimmed.fq", 'rb')
      records = SeqIO.parse(f_in, 'fastq')
      count=0
      for r in records:
        count+=1
        r.id= s.sampleID + "_" + str(count) # This is where I need the sample associated with it. 
        SeqIO.write(r, fna_out, 'fasta') # write it out to the file
      f_in.close()
    fna_out.close()





def print_help():
  print "\nUsage: trim_primers.py -- <input parameters>"
  print "Input order is unimportant. The '--' is required.\n"
  print "Mandatory inputs:"
  print "\t-i/--input-dir"
  print "\t\t An input directory of fastq files for trimming.\n"
  print "\t-o/--output-dir"
  print "\t\t An output directory that will contain the trimmed fastq files and combined .fna file\n"
  print "\t-f/--foward-primer"
  print "\t\t The sequence of the 5' primer typed inline. Do not specify a file.\n"
  print "\t-r/--reverse-primer"
  print "\t\t The sequence of the 3' primer typed inline. Do not specify a file.\n"
  print "\t-l/--length"
  print "\t\t The minimum length threshold for reads. Any reads shorter than this will be removed.\n"
  print "\t-q/--quality"
  print "\t\t The quality threshold for reads. All reads below this quality threshold will be removed.\n"


# creates the reverse compliment of a neucleotide sequence. 
# this is necessary for removing the 3' primer with cutadpt
def reverse_compliment(primer):
  # what are my characters and replacements?
  # C -> G
  # G -> C
  # A -> T
  # T -> A
  # Y -> R
  # R -> Y
  translation = {
    'c' : 'g',
    'C' : 'G',
    'g' : 'c',
    'G' : 'C',
    'a' : 't',
    'A' : 'T',
    't' : 'a',
    'T' : 'A',
    'y' : 'r',
    'Y' : 'R',
    'r' : 'y',
    'R' : 'Y'
  }
  output = ""
  print "input reverse primer:\n" + primer
  for i in range(len(primer)):
    output = translation.get(primer[i], "?") + output
  if output.find("?") != -1:
    raise Exception ("Could not interpret reverse primer sequence")
  print "output reverse compliment of reverse primer:\n" + output
  return output



if(__name__ == "__main__"):
  Trimmer()




