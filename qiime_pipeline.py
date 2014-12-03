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



class QiimePipeline(object):
  input = None
  output = None
  map = None
  reference = None
  num_cpus = 2
  min_seq_length = 200
  max_seq_length = 1000
  min_qual_score = 20
  qual_score_window = 50
  min_overlap = 20
  max_overlap = 100
  help = False




  def main(self):
    os.system("source /mnt/software/qiime/activate.sh")
    try:
      self.get_params()
      self.preprocess()
      print "The preprocessed.fna files are combined"
      os.system("cat " + self.output + "/preprocessed_fna/*trimmed_seqs.fna > " + self.output + "/preprocessed_fna/combined.fna")
      self.QiimeReport()
      
    except Exception as e:
      print e
      self.print_help()
      sys.exit()





  def get_params(self):
    letters = 'i:,o:,m:,r:,n:,l:,L:,q:,w:,p:,P:,h'
    keywords = ['input-dir=', 'output-dir=', 'map=', 'reference=', 'num_cpus=', 'min_seq_length=', 'max_seq_length=', 'min_qual_score=', 'qual_score_window=', 'min_overlap', 'max_overlap', 'help']
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
        else:
          raise Exception("Output directory already exists. Please specify a unique one.")
      elif o in ['-m', '--map']:
        if not os.path.exists(os.path.abspath(p)):
          raise Exception("Map file does not exit")
        self.map = os.path.abspath(p)
      elif o in ['-r', '--reference']:
	if not os.path.exists(os.path.abspath(p)):
	  raise Exception("Reference file does not exist.")
	self.reference = os.path.abspath(p)
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
      elif o in ['-q', '--min_qual_score']:
        if p.isdigit():
          self.min_qual_score = p
        else:
          raise Exception("Minimum quality score is not a number")
      elif o in ['-w', '--qual_score_window']:
        if p.isdigit():
          self.qual_score_window = p
        else:
          raise Exception("The quality score window is not a number")
      elif o in ['-p', '--min_overlap']:
        if p.isdigit():
          self.min_overlap = p
        else:
          raise Exception("Minimum overlap is not a number")
      elif o in ['-P', '--max_overlap']:
        if p.isdigit():
          self.max_overlap = p
        else:
          raise Exception("Maximum overlap is not a number")      
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





  def preprocess(self):
    #read mapping file
    f1 = open(self.map)
    lines=f1.readlines()
    n=2 
    while n < len(lines):
      self.sampleid = lines[n].split('\t')[0]
      self.fprimer = lines[n].split('\t')[2]
      self.treatment = lines[n].split('\t')[3]
      self.filename = lines[n].split('\t')[4]
      self.fileid = self.filename.split('.')[0]
      fileformat = self.filename.split('.')[-1]
      fileformat2 = self.filename.split('.')[-2]
      self.pairdreads = lines[n].split('\t')[5]
      self.rprimer = lines[n].split('\t')[6]
      self.fprimer_rc = lines[n].split('\t')[7]
      self.rprimer_rc = lines[n].split('\t')[8]
      self.readnumber = lines[n].split('\t')[9]
        
      #split mapping file
      if not os.path.exists(self.output + "/split_maps/"):
        os.makedirs(self.output + "/split_maps/")
      self.singlemap = self.sampleid + '_map.txt'         
      f2 = open (self.output + "/split_maps/" + self.singlemap, 'w')
      f2.write(lines[0])                       
      f2.write(lines[n])    
      f2.close()

      if not os.path.exists(self.output + "/" + self.sampleid + "/converted_seqs/"):
        os.makedirs(self.output + "/" + self.sampleid + "/converted_seqs/")
      if not os.path.exists(self.output + "/preprocessed_fna/"):
        os.makedirs(self.output + "/preprocessed_fna/")
      if not os.path.exists(self.output + "/preprocessing_logs"):
        os.makedirs(self.output + "/preprocessing_logs") 
      
      if fileformat == 'sff':
        print (str(n-1)+ " " + self.filename + ' is being denoised...')
        self.denoise_454()

      elif fileformat == 'bam':
        print (str(n-1)+ " " + self.filename + ' is being converted...')
        self.filter_pgm()
	
      elif fileformat == 'fastq' or fileformat == 'fq':
        print (str(n-1)+ " " + self.filename + ' is being split...')
        SeqIO.convert(self.input + "/" + self.filename, "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.qual", "qual")
        SeqIO.convert(self.input + "/" + self.filename, "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.fasta", "fasta")
        self.filter()
        
      elif fileformat2 == 'fastq' and fileformat == 'gz' and self.readnumber == 'all':
        print (str(n-1)+ " Merging " + self.filename + " and " + self.pairdreads)
        self.filter_miseq()

      elif fileformat2 == 'fastq' and fileformat == 'gz':
        print (str(n-1)+ " Merging " + self.filename + " and " + self.pairdreads)
        self.filter_miseq_group()

      else:
        print 'Input file format is invalid. Please check your mapping file.'
      n+=1
    f1.close()
    print "Data preprocessing is finished."  








  def denoise_454(self):
      
    print ("  Processing the .sff file...")
    os.system("process_sff.py -f -i " + self.input + "/" + self.filename + " -o " + self.output + "/" + self.sampleid + "/converted_seqs")

    print "  Running split_libraries..." 
    split_lib_command = "split_libraries.py -b 0 -p -g -o " + self.output + "/" + self.sampleid + "/slout -f " + self.output + "/" + self.sampleid + "/converted_seqs/" + self.fileid + ".fna -q " + self.output + "/" + self.sampleid + "/converted_seqs/" + self.fileid + ".qual -m " + self.output + "/split_maps/" + self.singlemap + " -l " + str(self.min_seq_length) + " -L " + str(self.max_seq_length) + " -s " + str(self.min_qual_score) + " -w " + str(self.qual_score_window)
    os.system(split_lib_command)
    os.system("cp " + self.output + "/" + self.sampleid + "/slout/split_library_log.txt " + self.output + "/preprocessing_logs/" + self.sampleid + "_slout_log")
    
    print "  Running denoise_wrapper..."
    denoise_wrapper_command = "denoise_wrapper.py -v -i " + self.output + "/" + self.sampleid + "/converted_seqs/" + self.fileid + ".txt -f " + self.output + "/" + self.sampleid + "/slout/seqs.fna -o " + self.output + "/" + self.sampleid + "/dwout/ -m " + self.output + "/split_maps/" + self.singlemap + " --force_overwrite"
    os.system(denoise_wrapper_command)

    print "  Running inflate_denoiser..."
    inflate_denoiser_command = "inflate_denoiser_output.py -c " + self.output + "/" + self.sampleid + "/dwout/centroids.fasta -s " + self.output + "/" + self.sampleid + "/dwout/singletons.fasta -f " + self.output + "/" + self.sampleid + "/slout/seqs.fna -d " + self.output + "/" + self.sampleid + "/dwout/denoiser_mapping.txt -o " + self.output + '/' + self.sampleid + "/denoiser.fna"
    os.system(inflate_denoiser_command)

    print "  Trimming primer sequences..."
    if not os.path.exists(self.output + "/" + self.sampleid + "/trimmed_seqs/"):
      os.makedirs(self.output + "/" + self.sampleid + "/trimmed_seqs/")
    log = self.output + "/" + self.sampleid + "/trimmed_seqs/cutadapt_log"
    (open(log, 'w')).close()
    os.system("cutadapt -b " +  self.fprimer + " -o " + self.output + "/" + self.sampleid + "/trimmed_seqs/ftrimmed.fna " + self.output + "/" + self.sampleid + "/denoiser.fna > " + log)
    os.system("cutadapt -b " +  self.rprimer + " -o " + self.output + "/" + self.sampleid + "/trimmed_seqs/rtrimmed.fna " + self.output + "/" + self.sampleid + "/trimmed_seqs/ftrimmed.fna >> " + log)
    os.system("cutadapt -b " +  self.fprimer_rc + " -o " + self.output + "/" + self.sampleid + "/trimmed_seqs/frctrimmed.fna " + self.output + "/" + self.sampleid + "/trimmed_seqs/rtrimmed.fna >> " + log)
    os.system("cutadapt -m " + str(self.min_seq_length) + " -b " +  self.rprimer_rc + " -o " + self.output + "/preprocessed_fna/" + self.sampleid + "_trimmed_seqs.fna " + self.output + "/" + self.sampleid + "/trimmed_seqs/frctrimmed.fna >> " + log)
    os.system("cp " + log + " " + self.output + "/preprocessing_logs/" + self.sampleid + "_cutadapt_log")         




  def filter(self):
	  
    print "  Running split_libraries..." 
    split_lib_command = "split_libraries.py -b 0 -p -g -o " + self.output + "/" + self.sampleid + "/slout -f " + self.output + "/" + self.sampleid + "/converted_seqs/seqs.fasta -q " + self.output + "/" + self.sampleid + "/converted_seqs/seqs.qual -m " + self.output + "/split_maps/" + self.singlemap + " -l " + str(self.min_seq_length) + " -L " + str(self.max_seq_length) + " -s " + str(self.min_qual_score) + " -w " + str(self.qual_score_window)
    os.system(split_lib_command)
    os.system("cp " + self.output + "/" + self.sampleid + "/slout/split_library_log.txt " + self.output + "/preprocessing_logs/" + self.sampleid + "_slout_log")
    os.system("cp " + self.output + "/" + self.sampleid + "/slout/histograms.txt " + self.output + "/preprocessing_logs/" + self.sampleid + "_histograms.txt")
        
    print "  Trimming primer sequences..."
    if not os.path.exists(self.output + "/" + self.sampleid + "/trimmed_seqs/"):
      os.makedirs(self.output + "/" + self.sampleid + "/trimmed_seqs/")
    log = self.output + "/" + self.sampleid + "/trimmed_seqs/cutadapt_log"
    (open(log, 'w')).close()
    os.system("cutadapt -b " +  self.fprimer + " -o " + self.output + "/" + self.sampleid + "/trimmed_seqs/ftrimmed.fna " + self.output + "/" + self.sampleid + "/slout/seqs.fna > " + log)
    os.system("cutadapt -b " +  self.rprimer + " -o " + self.output + "/" + self.sampleid + "/trimmed_seqs/rtrimmed.fna " + self.output + "/" + self.sampleid + "/trimmed_seqs/ftrimmed.fna >> " + log)
    os.system("cutadapt -b " +  self.fprimer_rc + " -o " + self.output + "/" + self.sampleid + "/trimmed_seqs/frctrimmed.fna " + self.output + "/" + self.sampleid + "/trimmed_seqs/rtrimmed.fna >> " + log)
    os.system("cutadapt -m " + str(self.min_seq_length) + " -b " +  self.rprimer_rc + " -o " + self.output + "/preprocessed_fna/" + self.sampleid + "_trimmed_seqs.fna " + self.output + "/" + self.sampleid + "/trimmed_seqs/frctrimmed.fna >> " + log)
    os.system("cp " + log + " " + self.output + "/preprocessing_logs/" + self.sampleid + "_cutadapt_log")






  def filter_pgm(self):

    os.system("/mnt/software/tophat2/bam2fastx --fastq " + self.input + "/" + self.filename + " >> " + self.output + "/" + self.sampleid + "/converted_seqs/seqs.fastq")
    SeqIO.convert(self.output + "/" + self.sampleid + "/converted_seqs/seqs.fastq", "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.qual", "qual")
    SeqIO.convert(self.output + "/" + self.sampleid + "/converted_seqs/seqs.fastq", "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.fasta", "fasta")
    self.filter()





  def filter_miseq(self):

    if not os.path.exists(self.output + "/" + self.sampleid + "/merged_fastq/"):
      os.makedirs(self.output + "/" + self.sampleid + "/merged_fastq/")
    flash_log = self.output + "/" + self.sampleid + "/merged_fastq/flash_log"
    (open(flash_log, 'w')).close()
    os.system ("/mnt/software/flash " + self.input + "/" + self.filename + " " + self.input + "/" + self.pairdreads + " -d " + self.output + "/" + self.sampleid + "/merged_fastq/ -o flash -m " + str(self.min_overlap) + " -M " + str(self.max_overlap) + " > " + flash_log)
    os.system("cp " + self.output + "/" + self.sampleid + "/merged_fastq/flash_log " + self.output + "/preprocessing_logs/" + self.sampleid + "_flash_log")
    
    SeqIO.convert(self.output + "/" + self.sampleid + "/merged_fastq/flash.extendedFrags.fastq", "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.qual", "qual")
    SeqIO.convert(self.output + "/" + self.sampleid + "/merged_fastq/flash.extendedFrags.fastq", "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.fasta", "fasta")

    self.filter()






  def filter_miseq_group(self):

    if not os.path.exists(self.output + "/" + self.sampleid + "/merged_fastq/"):
      os.makedirs(self.output + "/" + self.sampleid + "/merged_fastq/")
    flash_log = self.output + "/" + self.sampleid + "/merged_fastq/flash_log"
    (open(flash_log, 'w')).close()
    os.system ("/mnt/software/flash " + self.input + "/" + self.filename + " " + self.input + "/" + self.pairdreads + " -d " + self.output + "/" + self.sampleid + "/merged_fastq/ -o flash -m " + str(self.min_overlap) + " -M " + str(self.max_overlap) + " > " + flash_log)
    os.system("cp " + self.output + "/" + self.sampleid + "/merged_fastq/flash_log " + self.output + "/preprocessing_logs/" + self.sampleid + "_flash_log") 
    
    print "  Getting " + self.readnumber + " reads..."
    merged_fq = self.output + "/" + self.sampleid + "/merged_fastq/flash.extendedFrags.fastq"
    group_fq = self.output + "/" + self.sampleid + "/merged_fastq/group_" + self.readnumber + ".fastq"
    f1=open(merged_fq,"r")
    f2= open(group_fq,"w")
    lines=f1.readlines()
    for i in range (0, 4*int(self.readnumber)):
      f2.write(lines[i])
    f2.close()

    SeqIO.convert(group_fq, "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.qual", "qual")
    SeqIO.convert(group_fq, "fastq", self.output + "/" + self.sampleid + "/converted_seqs/seqs.fasta", "fasta")

    self.filter()	






  def QiimeReport(self):
	  
    QiimeReport_command = "python " + os.path.join(script_path, "qiime_report.py") + " -i " + self.output + "/preprocessed_fna/combined.fna -o " + self.output + " -m " + self.map + " -n " + str(self.num_cpus)
    if (not (self.reference is None)):
      QiimeReport_command += " -r " + self.reference
    os.system(QiimeReport_command)
         
 



  def __init__(self):
    self.main()





  def print_help(self):
    print '''
    qiime_pipeline.py runs the qiime pipeline for 16S rRNA metagenomic analysis. It takes demultiplexed SFF, BAM, and FASTQ files from 454, pgm and illumina sequencers.
    Required Parameters:
    (-i, --input-dir) The input directory containing all the data you want to analyzed. 
    (-o, --output-dir) The output directory where you want the results to be located.
    (-m, --map) The tab-delimited qiime mapping file that contains metadata for the samples to be analyzed.
    Optional Parameters:
    (-r, --reference) Reference database for OTU clustering. It should be in fasta format. If a reference is given, pick_open_reference_otus.py runs. If not, pick_de_novo_otus.py runs.
    (-n, --num_cpus) The number of cores to use for the parallel portions of the pipeline. [default: 2]
    (-l, --min-seq-length) Minimum sequence length, in nucleotides. [default: 200]
    (-L, --max-seq-length) Maximum sequence length, in nucleotides. [default: 1000]
    (-q, --min-qual-score) Min average qual score allowed in reads. [default: 20]
    (-w, --qual_score_window) Enable sliding window test of quality scores. If the average score of a continuous set of w nucleotides falls below the threshold, the sequence is discarded.[default: 50]
    (-p, --min_overlap) The minimum required overlap length between two reads to provide a confident overlap. This option is used for merging miseq paied-end reads. [Default: 20bp]
    (-P, --max_overlap) Maximum overlap length expected in approximately 90% of read pairs. This option is used for merging miseq paied-end reads. [Default: 100bp]
    Other options: 
    (-h, --help) Display this help dialogue and exit.
    For complete information about this pipeline, please refer to the manual. 
    '''




if(__name__ == "__main__"):
  QiimePipeline()
