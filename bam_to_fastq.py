#!/mnt/software/epd/bin/ipython
import os
import getopt
import sys
import run_qiime # I'm going to comment this out and rerun it. 
import collections
import csv

class Converter(object):
  input = None
  output = None
  map = None
  input_files = None
  #output_files = None
  samples = []

  def get_params(self):
    letters = 'i:,o:,m:'
    keywords = ['input-dir=', 'output-dir=', 'map=']

    opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)
    for o,p in opts:
      if o in ['-i', '--input-dir']:
        self.input = p
      elif o in ['-o', '--output-dir']:
        self.output = p
      elif o in ['-m', '--map']:
        self.map = p
      else:
        raise Exception("Improper Input to converter")
    if(self.input is None or self.output is None or self.map is None):
      raise Eception("Incomplete Input") 


  # What I'm about to do here:
  # I need to modify this so that it looks through the map file and gets the input file names
  # it will then relate these input filenames to the sample IDs and output the converted files into their new names based on sample id
  def convert_to_fastq(self):
    samples = parse_map(self.map)
    os.system("source /mnt/software/qiime/qiime_install/activate.sh")
    if(not os.path.exists(self.input)):
      raise Exception("Invalid Input Path")
   # check that the destination directory exists, if not, make it
    if (not os.path.exists(self.output + "/fastq/")):
      os.makedirs(self.output + "/fastq/")
      print "Directory " + self.output + "does not exist. Creating it now."
 
    # iterate through bam files calling bam2fastx
    # I can parallelize this for loop to speed this up. 
    for s in samples:
      print "\nConverting " + s.filename + " to " + s.sampleID + ".fq"
      os.system("bam2fastx --fastq " + str(self.input) + "/" + s.filename + " >> " + str(self.output) + "/fastq/" + s.sampleID + ".fq")
      print "Done converting " + s.sampleID

  def __init__(self):
    self.main()

  def main(self):
    self.get_params()
    self.samples = run_qiime.parse_map(self.map)
    self.convert_to_fastq()
    print "\nConversion Complete"



  # This should look through the map file and pull out the following information:
  # forward and reverse primer sequences. 
  # input file names
  # group ID
  #how to do this:
  # read first line and associate indexes with values. so sample ID = 0 etc
  #loop through all the rows, pulling out data. I might want to store everything as records where each sample gets a record associated with it and we store these in a sample record list. 
  # the record should have an additional parameters list that sores everything I don't have a specific value for. we might need them later and it's just nice to have. 
  # after this, I'm going to need to go through and edit a bunch of the code elsewhere. 
def parse_map(map):
  fprimer_index = 2
  rprimer_index = -1
  group_index = -1
  filename_index = -1
  sampleID_index = 0
  barcode_index = 1
  description_index = -1
  other_indexes = []
  sample_list =[]
  Sample = collections.namedtuple('Sample', ['fprimer', 'rprimer', 'group', 'filename', 'sampleID', 'barcode', 'description', 'other'])
  with open(map, "rb") as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    first_line = reader.next()
    for i in range(len(first_line)):
      if first_line[i].lower() == "groupid":
        group_index = i
      elif first_line[i].lower() == "reverseprimer":
        rprimer_index = i
      elif first_line[i].lower() == "filename":
        filename_index = i
      elif first_line[i].lower() == "description":
        description_index = i
      else:
        if i!=fprimer_index and i!=sampleID_index and i!=barcode_index:
          other_indexes += [i]

    for row in reader:
      if row[0][0] != '#':
        sample_data = Sample(fprimer = row[fprimer_index], rprimer = row[rprimer_index], group = row[group_index], filename = row[filename_index],
          sampleID = row[sampleID_index], barcode = row[barcode_index], description = row[description_index], other = [])
        for i in other_indexes:
          sample_data = sample_data._replace(other= sample_data.other + [row[i]])
        #print str(sample_data);
        sample_list += [sample_data]
  return sample_list

if(__name__ == "__main__"):
  Converter()



