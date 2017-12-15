#!/usr/bin/python
##usage: python modify_snpEff_VCF.py input_variants.vcf.gz output_file.vcf
##this script modifies the ANN,LOF and NMD fields from snpEff VCF to resemble a more standard VCF 
##that can be processed through bior_vcf2xls.pl

import gzip
import sys

f=open(sys.argv[2],'w')
counter=0
longestanno=[]
if str(sys.argv[1]).endswith('.gz'):
	vcf=gzip.open(sys.argv[1],'r')
else:
	vcf=open(sys.argv[1],'r')
	for line in vcf:
		if line.startswith('#'):
			f.write(line)
		else:
			part1='\t'.join(line.rstrip().split('\t')[0:7])
			lastpart='\t'.join(line.rstrip().split('\t')[8:12])
			##get INFO fields
			INFO=line.rstrip().split('\t')[7].split(';')
			#print INFO
			annoFields=["Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA_coord","CDS_coord","AA_coord","Distance","notes"]
			lof_nmd=["Gene_Name","Gene_ID","Number_of_transcripts_in_gene","Percent_of_transcripts_affected"]
			partINFO=''
			anno=''
			lof=''
			nmd=''
			for i in range(0,len(INFO)):
				#if (not INFO[i].startswith('ANN=')):
				#	partINFO+=INFO[i]+";"
				if (INFO[i].startswith('ANN=')):
					ann=INFO[i].rstrip().split('=')[1].split(',')
					for b in range(0,len(ann)):
						a=0
						while (a<16):
							anno+=str(annoFields[a])+"_"+str(b)+"="
							if (ann[b].rstrip().split('|')[a]):
								anno+=str(ann[b].rstrip().split('|')[a])+";"
							else:
								anno+="."+";"
							a+=1
					if (len(anno.split(';'))>counter):
						counter=len(anno.split(';'))
						longestanno=anno.split(';')
					#print anno
				elif (INFO[i].startswith('LOF=')):
					#print "blah"
					lofvals=INFO[i].rstrip().split('=')[1].replace("(","").replace(")","").split('|')
					#print lofvals
					for c in range(0,4):
						lof+=str(lof_nmd[c])+"_LOF="
						if (lofvals[c]):
							lof+=str(lofvals[c])+";"
						else:
							lof+=".;"
				elif (INFO[i].startswith('NMD=')):
					nmdvals=INFO[i].rstrip().split('=')[1].replace("(","").replace(")","").split('|')
					for d in range(0,4):
						nmd+=str(lof_nmd[d])+"_NMD="
						if (nmdvals[d]):
							nmd+=str(nmdvals[d])+";"
						else:
							nmd+=".;"
				else:
					partINFO+=INFO[i]+";"
			f.write(part1+"\t"+partINFO+(anno+lof+nmd).rstrip(';')+"\t"+lastpart+"\n")

f.close()
print "The highest number of annotations for any variant are"+str(counter) 
print ' '.join(longestanno)
