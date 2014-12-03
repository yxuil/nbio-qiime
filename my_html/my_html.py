import os
import time
import trim_primers
import glob

outputDir= "../test_output/"
map = "../../../steele_data/steele_map_corrected.txt"
PIName = "PI Name"
ref_db = "greengenes"


os.system("mkdir -p "+outputDir+"/report_files/images")
f=open(outputDir+"/report.html","w")
f.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n")
f.write("<html>\n")
f.write("<head>\n")
f.write("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n")

f.write("<style type=\"text/css\">\n")
f.write("body, th, td, p { font-family: helvetica, sans-serif; font-size: small }\n")
f.write("h1 { font-size: medium; text-align: center }\n")
f.write("h2, h3 { font-size: small; text-align: left }\n")
f.write("table {}\n") # I could change this to get rid of the table border. make two types of table, one with a border and the default
# without. then add tags to all the tables I use so that they have the borders and then the qiime tables will follow the default table behavior. 
f.write(".brc { border: 1px solid black; border-collapse: collapse; margin: 10px 10px 15px 0; width: 820px }")

#f.write("table tr:nth-child(even) td { background-color: #eeeeee }\n")# uncommenting this fixes it.
f.write("th { text-align: left; color: white; background-color: #aa3333; padding: 3px 15px 3px 3px; border-bottom: 1px solid black; width: 150px }\n")
f.write("td { padding: 3px 3px 3px 4px; border-bottom: 1px solid gray  }\n")
f.write(".ntitle { color: black; font-family:Arial,Verdana; font-size:11; font-weight:bold;}")
f.write(".left { text-align: left }\n")
f.write(".right { text-align: right }\n")
f.write(".center { text-align: center }\n")
f.write(".boldcentered { text-align: center; font-weight: bold }\n")
f.write("</style>\n")
#f.write("<link rel=\"stylesheet\" href=\"./css/qiime_style.css\" type=\"text/css\">")




f.write("</head>\n")
f.write("<body>\n")
f.write("<div class=\"boldcentered\">UW Biotechnology Center Bioinformatics Resource Center - "+time.asctime()+"</div>\n")
f.write('<div style="margin-top: 10px" class="center">Please direct questions to: <strong>brc@biotech.wisc.edu</strong></div>'+"\n");
f.write("<hr>\n");
f.write("<ol>\n");
f.write("<li><a href='#summary'>Summary</a></li>\n");
f.write("<li><a href='#workflow'>Work Flow</a></li>\n");
f.write("<li><a href='#taxa_summary'>Taxa Summary</a></li>\n");
f.write("<li><a href='#arare'>Alpha Rarefaction</a></li>\n");
f.write("<li><a href='#bdiv'>Beta Diversity</a></li>\n");
f.write("<li><a href='#phylogeny'>Phylogenetic Tree</a></li>\n")
f.write("<li><a href='#heatmap'>OTU Heatmap</a></li>\n")
f.write("<li><a href='#software'>Software</a></li>\n");
f.write("</ol>\n");
f.write("<hr>\n");



samples = trim_primers.parse_map(map)
#summary
f.write("<div class=\"left\"><strong><a id=summary name=summary>Summary</a></strong></div>")
f.write("<table class=\"brc\">\n")
#f.write("<caption>Summary</caption>\n")
f.write("<tr><th scope='row'>PI Name</th><td>"+PIName+"</td></tr>\n")
f.write("<tr><th scope='row'>Project Description</th><td>"+samples[0].description+"</td></tr>\n")
sample_names = ""
for s in samples:
  sample_names += s.sampleID + ", "
f.write("<tr><th scope='row'>Samples</th><td>"+sample_names+"</td></tr>\n")
if(not (ref_db is None)):
  f.write("<tr><th scope='row'>Reference Database</th><td>"+ref_db+"</td></tr>\n")
#f.write("<tr><th>Report Generation Date</th><td>"+time.asctime()+"</td></tr>\n")
f.write("</table>\n")
f.write("<hr>\n");
#reads and coverage
'''
f.write("<div class='left'><strong><a id='reads' name='reads'>Reads statistics</a></strong> (click the links in mapped reads field to download the final bam files)</div>")
f.write("<table>\n")
f.write("<tr>\n")
f.write("<th scope='col'>Sample Name</th>\n")
f.write("<th scope='col'>Total Reads</th>\n")
f.write("<th scope='col'>Mapped Reads</th>\n")
f.write("<th scope='col'>Mapped Reads Percentage</th>\n")
f.write("<th scope='col'>Coverage</th>\n")
'''
# So I need to figure out where to find the read statistics. 


# The workflow
f.write('<hr>\n')
f.write("<div class='left'><strong><a id='workflow' name='workflow'>Work flow chart</a></strong></div><table>\n")
f.write('<tr>\n')
f.write('<td class=center><img src="../my_html/Workflow.png" alt="Work flow chart"></td></tr>\n')
f.write('</table>\n')
f.write('<hr>\n')




# The taxonomy summary
# How about I put a table of all the stuff in the excell sheet and include the pictures. Basically just replicate the html page that the pipeline creates. 
# I could just link to it but it would look much worse, let's see how much work it'd be to copy over the stuff from the html that qiime generates. 
# okay. I think I can just hard code the header and scripts from the qiime output then just grab everything after the line containing "otu_table_L5.txt"

# how do I make this happen? 
# I'll have to use a filereader and iterate though line by line until I find the string I'm looking for
# or there might be a find method that returns line numbers. 

path = os.path.abspath(outputDir + 'otus/taxa_summary/taxa_summary_plots/bar_charts.html')
f.write("<div class=\"left\"><strong><a id=taxa_summary name=taxa_summary>Taxa Summary</a></strong></div>")
f.write("<li><a href=\'" + path + "\'>View Full Interactive Taxa Summary</a></li>\n");

with open(outputDir + 'otus/taxa_summary/taxa_summary_plots/bar_charts.html', 'r') as bar_summary:
  start=False
  for line in bar_summary:
    if start:
      str = line.replace('a href=\'', 'a href=\'' + 'otus/taxa_summary/taxa_summary_plots/')
      str = str.replace('a href=\"', 'a href=\"' + 'otus/taxa_summary/taxa_summary_plots/')
      str = str.replace('img src=\'', 'img src=\'' + 'otus/taxa_summary/taxa_summary_plots/')
      f.write(str)
    if line.find("otu_table_L5.txt") != -1:
      start = True

f.write("<hr>\n");


f.write("<div class=\"left\"><strong><a id=arare name=arare>Alpha Rarefaction</a></strong></div>")
path = os.path.abspath(outputDir + 'arare/alpha_rarefaction_plots/rarefaction_plots.html')
f.write("<li><a href=\'" + path + "\'>View Alpha Rarefaction Plots</a></li>\n");

f.write("<hr>\n")
f.write("<div class=\"left\"><strong><a id=bdiv name=bdiv>Beta Diversity</a></strong></div>")


path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/unweighted_unifrac_2d_continuous/unweighted_unifrac_pc_2D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 2D Continuous Unweighted Unifrac </a></li>\n");

path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/unweighted_unifrac_2d_discrete/unweighted_unifrac_pc_2D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 2D Discrete Unweighted Unifrac </a></li>\n");

path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/weighted_unifrac_2d_continuous/weighted_unifrac_pc_2D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 2D Continuous Weighted Unifrac </a></li>\n");

path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/weighted_unifrac_2d_discrete/weighted_unifrac_pc_2D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 2D Discrete Weighted Unifrac </a></li>\n");


path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/unweighted_unifrac_3d_continuous/unweighted_unifrac_pc_3D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 3D Continuous Unweighted Unifrac</a></li>\n");

path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/unweighted_unifrac_3d_discrete/unweighted_unifrac_pc_3D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 3D Discrete Unweighted Unifrac </a></li>\n");

path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/weighted_unifrac_3d_continuous/weighted_unifrac_pc_3D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 3D Continuous Weighted Unifrac </a></li>\n");

path = os.path.abspath(glob.glob(outputDir + '*bdiv*')[0] + '/weighted_unifrac_3d_discrete/weighted_unifrac_pc_3D_PCoA_plots.html')
f.write("<li><a href=\'" + path + "\'>View 3D Discrete Weighted Unifrac </a></li>\n");

f.write("<hr>\n");


#Heatmap
f.write("<div class='left'><strong><a id='heatmap' name='heatmap'>Heatmap</a></strong></div>")
path = os.path.abspath(outputDir + 'otus/heatmap/otu_table.html')
f.write("<li><a href=\'" + path + "\'>View OTU Heatmap</a></li>\n");
f.write("<hr>\n")


# Phylogeny
f.write("<div class='left'><strong><a id='phylogeny' name='phylogeny'>Phylogenetic Tree</a></strong></div>")
path = os.path.abspath(outputDir + 'otus/rep_set.tre')
f.write("<li><a href=\'" + path + "\'>Right click and select \"save as\" to download tree</a></li>\n");
f.write("<hr>\n")



#software
f.write("<div class='left'><strong><a id='software' name='software'>Software</a></strong></div>")
f.write("<table class=\"brc\">\n")
f.write("<tr><th>Main Pipeline</th><td>QIIME</td></tr>\n")
f.write("<tr><th>OTU Clustering</th><td>Uclust</td></tr>\n")
f.write("<tr><th>Quality Filtering</th><td>Usearch</td></tr>\n")
f.write("<tr><th>Alignment</th><td>PyNAST</td></tr>\n")
f.write("<tr><th>Taxonmy Assignment</th><td>RDP Classifier</td></tr>\n")
f.write("<tr><th>Phylogeny Generation</th><td>FastTree</td></tr>\n")
f.write("</table>\n")
f.write("<hr>\n");





'''
if(len(TargetRegion)>0):
    f.write("<th scope='col'>Coverage in Target Region</th>\n")
f.write("</tr>\n")
for i in range(len(Samples)):
    DeliveryDir=configs[i].get('Sample',"DeliveryDir")
    SampleName=configs[i].get('Sample','SampleName')
    f.write("<tr>\n")
    f.write("<td>"+configs[i].get('Sample','SampleName')+"</td>\n")
    f.write("<td>"+configs[i].get('Sample','TotalReads')+"</td>\n")
    f.write("<td><a href='./"+SampleName+"/chrAll.bam'>"+configs[i].get('Sample','MappedPassed')+"</a></td>\n")
    if(int(configs[i].get('Sample','TotalReads'))==0):
        f.write("<td>0.0%</td>\n")
    else:
#           f.write("<td>"+"%3.1f" % (float(configs[i].get('Sample','MappedPassed'))/int(configs[i].get('Sample','TotalReads'))*100)+"%</td>\n")
        f.write("<td>"+"%3.1f" % (float(configs[i].get('Sample','MappedPassed'))/int(configs[i].get('Sample','TotalReadsPassed'))*100)+"%</td>\n")
    f.write("<td><a href='#dcov'>"+configs[i].get('Sample','AverageCoverageAll')+"</a></td>\n")
    if(len(TargetRegion)>0):
        f.write("<td><a href='#dcov'>"+configs[i].get('Sample','AverageCoverageTargeted')+"</a></td>\n")
    f.write("</tr>\n")
f.write("</table>\n")
f.write("<hr>\n");
'''
'''
    #variants
    f.write("<div class='left'><strong><a id=variants name=variants>Variants statistics</a></strong> (click the links in the table to download variants files)</div>")
    f.write("<div class='left'><strong>GATK calls</strong></div>")
    f.write("<table>\n")
    #f.write("<caption>Variants</caption>\n")
    f.write("<tr>\n")
    f.write("<th>Sample Name</th>\n")
    f.write("<th>Total Variants</th>\n")
    f.write("<th>SNPs</th>\n")
    f.write("<th>Indels</th>\n")
    if(len(TargetRegion)>0):
        f.write("<th>Total Variants in Target Region</th>\n")
        f.write("<th>SNPs in Target Region</th>\n")
        f.write("<th>Indels in Target Region</th>\n")
    f.write("</tr>\n")
    for i in range(len(Samples)):
        DeliveryDir=configs[i].get('Sample',"DeliveryDir")
        SampleName=configs[i].get('Sample','SampleName')
        VCFDir="./"+SampleName+"/vcf"
        allsnv=VCFDir+"/all.vcf"
        allsnp=VCFDir+"/all.snp.vcf"
        allindel=VCFDir+"/all.indel.vcf"
        alltargeted=VCFDir+"/all.targeted.vcf"
        alltargetedsnp=VCFDir+"/all.targeted.snp.vcf"
        alltargetedindel=VCFDir+"/all.targeted.indel.vcf"
        f.write("<tr>\n")
        f.write("<td>"+configs[i].get('Sample','SampleName')+"</td>\n")
        f.write("<td><a href='"+allsnv+"'>"+configs[i].get('Sample','VariantAll')+"</a></td>\n")
        f.write("<td><a href='"+allsnp+"'>"+configs[i].get('Sample','SNPAll')+"</a></td>\n")
        f.write("<td><a href='"+allindel+"'>"+configs[i].get('Sample','IndelAll')+"</a></td>\n")
        if(len(TargetRegion)>0):
            f.write("<td><a href='"+alltargeted+"'>"+configs[i].get('Sample','TargetedAll')+"</a></td>\n")
            f.write("<td><a href='"+alltargetedsnp+"'>"+configs[i].get('Sample','TargetedSNPAll')+"</a></td>\n")
            f.write("<td><a href='"+alltargetedindel+"'>"+configs[i].get('Sample','TargetedIndelAll')+"</a></td>\n")
        f.write("</tr>\n")
    f.write("</table>\n")
    f.write("<div class='left'><strong>Pindel calls</strong></div>")
    f.write("<table>\n")
    #f.write("<caption>Variants</caption>\n")
    f.write("<tr>\n")
    f.write("<th>Sample Name</th>\n")
    f.write("<th>Deletions</th>\n")
    f.write("<th>Small Insertions</th>\n")
    f.write("<th>Large Insertions</th>\n")
    f.write("<th>Break Points</th>\n")
    f.write("<th>Tandem Duplicatess</th>\n")
    f.write("<th>Inversions</th>\n")
    f.write("</tr>\n")
    for i in range(len(Samples)):
        ReportDir=configs[i].get('Sample',"ReportDir")
        SampleName=configs[i].get('Sample','SampleName')
        pindel_dir="./"+SampleName+"/pindel"
        chrAll_D=pindel_dir+"/chrAll_D"
        chrAll_D_vcf=pindel_dir+"/chrAll_D.vcf"
        chrAll_SI=pindel_dir+"/chrAll_SI"
        chrAll_SI_vcf=pindel_dir+"/chrAll_SI.vcf"
        chrAll_INV=pindel_dir+"/chrAll_INV"
        chrAll_INV_vcf=pindel_dir+"/chrAll_INV.vcf"
        chrAll_TD=pindel_dir+"/chrAll_TD"
        chrAll_TD_vcf=pindel_dir+"/chrAll_TD.vcf"
        chrAll_LI=pindel_dir+"/chrAll_LI"
        chrAll_LI_vcf=pindel_dir+"/chrAll_LI.vcf"
        chrAll_BP=pindel_dir+"/chrAll_BP"
        chrAll_BP_vcf=pindel_dir+"/chrAll_BP.vcf"
        f.write("<tr>\n")
        f.write("<td>"+configs[i].get('Sample','SampleName')+"</td>\n")
        f.write("<td><a href='"+chrAll_D+"'>"+configs[i].get('Sample','chrAll_D')+"</a></td>\n")
        f.write("<td><a href='"+chrAll_SI+"'>"+configs[i].get('Sample','chrAll_SI')+"</a></td>\n")
        f.write("<td><a href='"+chrAll_LI+"'>"+configs[i].get('Sample','chrAll_LI')+"</a></td>\n")
        f.write("<td><a href='"+chrAll_BP+"'>"+configs[i].get('Sample','chrAll_BP')+"</a></td>\n")
        f.write("<td><a href='"+chrAll_TD+"'>"+configs[i].get('Sample','chrAll_TD')+"</a></td>\n")
        f.write("<td><a href='"+chrAll_INV+"'>"+configs[i].get('Sample','chrAll_INV')+"</a></td>\n")
        f.write("</tr>\n")
    f.write("</table>\n")
    f.write("<hr>\n");
    #muTect calls
    if(config.has_section('muTect')):
        f.write("<div class='left'><strong><a id=mutect name=mutect>muTect calls statistics</a></strong> (click the links in the table to download called somatic mutation files, which can be opened by Excel)</div>")
        f.write("<table>\n")
        #f.write("<caption>muTect calls</caption>\n")
        f.write("<tr>\n")
        f.write("<th>Tumor Sample</th>\n")
        f.write("<th>Normal Sample</th>\n")
        f.write("<th>Total calls</th>\n")
        if(len(outintervalfile)>0):
            f.write("<th>Total calls in Target Region</th>\n")
        f.write("</tr>\n")
        for i in range(len(muTectWorkers)):
            DeliveryDir=muTectWorkers[i].DeliveryDir
            CancerSampleName=muTectWorkers[i].CancerConfig.get('Sample','SampleName')
            NormalSampleName=muTectWorkers[i].NormalConfig.get('Sample','SampleName')
            FileBase=muTectWorkers[i].FileBase
            muTectall="./muTect_"+FileBase+"/"+FileBase+".all.keep.csv"
            muTectint="./muTect_"+FileBase+"/"+FileBase+".int.keep.csv"
            f.write("<tr>\n")
            f.write("<td>"+CancerSampleName+"</td>\n")
            f.write("<td>"+NormalSampleName+"</td>\n")
            f.write("<td><a href='"+muTectall+"'>"+str(muTectWorkers[i].n_all)+"</a></td>\n")
            if(len(outintervalfile)>0):
                f.write("<td><a href='"+muTectint+"'>"+str(muTectWorkers[i].n_int)+"</a></td>\n")
            f.write("</tr>\n")
        f.write("</table>\n")
        f.write("<hr>\n");
    #software
    f.write("<div class='left'><strong><a id='software' name='software'>Software</a></strong></div>")
    f.write("<table>\n")
    #f.write("<caption>Software</caption>\n")
    #f.write("<th>Quality Control</th><td>FastX ToolKit</td>\n")
    f.write("<tr><th>Alignment</th><td>BWA</td></tr>\n")
    f.write("<tr><th>CleanUp</th><td>Picard & GATK</td></tr>\n")
    f.write("<tr><th>Variants Calling</th><td>GATK, Pindel</td></tr>\n")
    f.write("<tr><th>Sample Comparison</th><td>muTect</td></tr>\n")
    #f.write("<tr><th>Annotation</th><td>Annovar</td></tr>\n")
    f.write("</table>\n")
    f.write("<hr>\n");
    #databases
    f.write("<div class='left'><strong><a id='databases' name='databases'>Databases</a></strong></div>")
    f.write("<table>\n")
    #f.write("<caption>Databases</caption>\n")
    #f.write("<th>Quality Control</th><td>FastX ToolKit</td>\n")
    f.write("<tr><th>Reference Genome</th><td>"+config.get('Reference','GenomeName')+", "+config.get('Reference','BuildName')+"</td></tr>\n")
    SNP=config.get('Reference','SNP')
    if(len(SNP)>0):
        f.write("<tr><th>SNP</th><td>"+config.get('Reference','SNPDescription')+"</td></tr>\n")
    INDEL=config.get('Reference','Indel')
    if(len(INDEL)>0):
        f.write("<tr><th>Indel</th><td>"+config.get('Reference','IndelDescription')+"</td></tr>\n")
    #f.write("<tr><th>Annotation</th><td>"+config.get('Reference','Annovardb')+"</td></tr>\n")
    f.write("</table>\n")
    f.write("<hr>\n");
#work flow
    f.write("<div class='left'><strong><a id='workflow' name='workflow'>Work flow chart</a></strong></div>")
    f.write("<table>\n")
   #f.write("<caption>Flow Chart</caption>\n")
    f.write("<tr>\n")
    f.write("<td class=center>")
    os.system("cp "+config.get('Other','FlowChart')+" "+ReportDir+"/report_files/images/flowchart.png")
    f.write('<img src="./report_files/images/flowchart.png" alt="Work flow chart">')
    f.write("</td>")
    f.write("</tr>\n")
    f.write("</table>\n")
    f.write("<hr>\n");
#coverage chart
    dcov=config.get('Pipeline','dcov')
    cts=dcov.split(',')
    per=[]
    samplenames=[]
    f.write("<div class='left'><a id='dcov' name='dcov'><strong>Depth of coverage chart</strong></a></div>")
    f.write("<table>\n")
    #f.write("<caption>Coverage Chart</caption>\n")
    f.write("<tr>\n")
    f.write("<th>Coverage</th>\n")
    f.write(string.join(["<td>"+ct+"</td>" for ct in cts],"\n"))
    f.write("</tr>\n")
    cts=map(int,cts)
    for i in range(len(Samples)):
        if(configs[i].get('Sample','CoverageDataAvailableAll')=="1"):
            samplenameall=configs[i].get('Sample','SampleName')+":whole genome"
            samplenames.append(samplenameall)
            percentageall=configs[i].get('Sample','CoveragePercentageAll').strip().split(',')
            per.append(map(float,percentageall))
            f.write("<tr>\n")
            f.write("<th>"+samplenameall+"</th>\n")
            f.write(string.join(["<td>"+perc+"</td>" for perc in percentageall],"\n"))
            f.write("</tr>\n")
        if(len(TargetRegion)>0 and configs[i].get('Sample','CoverageDataAvailableTargeted')=="1"):
            samplenametargeted=configs[i].get('Sample','SampleName')+":target region"
            samplenames.append(samplenametargeted)
            percentagetargeted=configs[i].get('Sample','CoveragePercentageTargeted').strip().split(',')
            per.append(map(float,percentagetargeted))
            f.write("<tr>\n")
            f.write("<th>"+samplenametargeted+"</th>\n")
            f.write(string.join(["<td>"+perc+"</td>" for perc in percentagetargeted],"\n"))
            f.write("</tr>\n")
    f.write("</table>\n")
    DeliveryDir=configs[i].get('Sample',"DeliveryDir")
    pct_curve(per, samplenames, cts, ReportDir+"/report_files/images/dcov.png")
    f.write('<div class="left"><table><tr><td><img src="./report_files/images/dcov.png" alt="Depth of coverage chart"></td></tr></table></div>')
    '''
f.write("</body>\n")
f.write("</html>\n")
f.close()
print "Report done."


