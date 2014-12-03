#!/mnt/software/epd/bin/ipython

import os
import getopt
import sys
import csv
import collections
class QiimePipeline(object):
  fna = None
  cores = None
  map = None
  ref = None
  out = None
  sort = None
  # main class just calls the two methods to process the data
  def main(self):
    try:
      os.system("source /mnt/software/qiime/qiime_install/activate.sh")
      self.get_params()
      self.run_commands()
    except Exception as e:
      print e #uncommment for debugging or when running as a standalone app and not part of the pipeline.
      exit()

  # what do I need?
  # I need the location of the input fna file
  # I need an output directory
  # I need the # of cores on the machine (probably a built in function to figure that out. but I should probably prompt for the number they want to use at the beginning)
  # I need the location of the otu table
  # I need the location of the map
  # I need to know whether we're doing de novo or reference based analysis
  def get_params(self):
    letters = 'i:,o:,m:,r:,c:,s:' # letters for passing args
    keywords = ['input-fna=', 'output-dir=', 'map=', 'reference=', 'cores=', 'sort='] # keywords for passing args

    opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords) # start at second element because first is script name. opts, extraparams are lists of tuples
    for o,p in opts:
      if o in ['-i', '--input-fna']:
        self.fna = p
      elif o in ['-o', '--output-dir']:
        self.out = p
      elif o in ['-m', '--map']:
        self.map = p
      elif o in ['-r', '--reference']:
        self.ref = p
      elif o in ['-c', '--cores']:
        self.cores = p
      elif o in ['-s', '--sort']:
        self.sort = p
      else:
        raise Exception("Improper Input")

    #sanity check
    if(self.fna is None or self.out is None or self.map is None or self.sort is None):
      raise Exception("Imporper Input")
    if(self.cores is None):
      self.cores = 1
    if(not os.path.exists(self.fna)):
      raise Exception("Input Directory does not exist")



  def run_commands(self):
    # Generate OTUs
    print "Picking OTUs. This will take a long time."
    if(self.ref is None):
      #print ("pick_de_novo_otus.py -i " + self.fna + " -o " + self.out + "/otus -f -a -O " + str(self.cores) + ' -p pick_otu_params.txt')
      os.system("pick_de_novo_otus.py -i " + self.fna + " -o " + self.out + "/otus -f -a -O " + str(self.cores))# + ' -p pick_otu_params.txt')
    else:
      os.system("pick_open_reference_otus.py -i " + self.fna + " -o " + self.out + "/otus -r " + self.ref + " -f -a -O " + str(self.cores))

    # Check the OTU result
    print "Checking OTU result"
    os.system("print_biom_table_summary.py -i " + self.out + "/otus/otu_table.biom > " + self.out + "/otus/library_stats.txt")
    print "Finished creating OTUs"

    # sort the OTU table
    print "Sorting OTU table"
    os.system("sort_otu_table.py -i " + self.out + "/otus/otu_table.biom -s " + self.sort + " -m " + self.map + " -o " + self.out +"/otus/otu_table.biom")

    # make the heatmap
    print "Making OTU heatmap."
    os.system("make_otu_heatmap_html.py -i " + self.out + "/otus/otu_table.biom -o " + self.out + "/otus/heatmap")
    print "Finished creating heatmap"

    # make network for cytoscape
    print "Creating OTU network"
    os.system("make_otu_network.py -m " + self.map + " -i " + self.out + "/otus/otu_table.biom -o " + self.out + "/otus/network")
    print "Finished creating OTU network"

    # create taxa summary
    print "Creating taxa summary plots"
    os.system("summarize_taxa_through_plots.py -f -i " + self.out + "/otus/otu_table.biom -o " + self.out + "/otus/taxa_summary -m " + self.map)
    print "Finished creating taxa summary plots"

    #### In the future (like later today) I would like to make it so the user can specify his own parameters file.
    #### If the user doesn't specify a parameters file, it will default to ours
    # create alpha parameters
    print "Performing alpha diversity analysis and rarefation plots."
    if(not os.path.exists(self.out + "/arare/")):
      os.makedirs(self.out + "/arare/")
    f = open(self.out + "/arare/alpha_params.txt", 'w')
    f.write('')
    f.close()
    os.system("echo \"alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species\" > " + self.out + "/arare/alpha_params.txt")
    print "Finished creating alpha parameters"
    # create rarefaction curves
    os.system("alpha_rarefaction.py -f -i " + self.out + "/otus/otu_table.biom -m " + self.map + " -o " + self.out +
      "/arare/ -p " + self.out + "/arare/alpha_params.txt -t " + self.out + "/otus/rep_set.tre -a -O " +str(self.cores))
    print "finished generating rarefaction curves"

    # create beta diversity plots
    print "Creating beta diversity plots"   
    os.system("beta_diversity_through_plots.py -f -i " + self.out + "/otus/otu_table.biom -m " + self.map + " -o " + self.out + "/bdiv -t " + self.out + "/otus/rep_set.tre -e 18000")
    print "Finished creating beta diversity plots"


  def __init__(self):
    self.main()



# scan through the map file and create 'sample' objects that contain all the metadata about each sample. 
def parse_map(map):

  # Find the indexes of all metadata that doesn't have a definite location. 
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
    # read header line and figure out the indexes of the metadata
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

    # go through the rows of the map and generate Sample objects for each sample
    for row in reader:
      if row[0][0] != '#':
        sample_data = Sample(fprimer = row[fprimer_index], rprimer = row[rprimer_index], group = row[group_index], filename = row[filename_index],
          sampleID = row[sampleID_index], barcode = row[barcode_index], description = row[description_index], other = []) 
        for i in other_indexes:
          sample_data = sample_data._replace(other= sample_data.other + [row[i]])
        sample_list += [sample_data]
  # return a list containing the metadata for all the samples. 
  return sample_list







if(__name__ == "__main__"):
  QiimePipeline()





