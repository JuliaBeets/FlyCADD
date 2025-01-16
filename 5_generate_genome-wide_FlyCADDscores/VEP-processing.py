#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 08-08-18

It takes an VEP output file as input and returns a tab delimited, encoded file.
And adds additional annotations.
The annotated variants file contains annotations for:
	the chromosome, position, reference base, alternative base, whenever it is a transversion, 
	GC%, CpG%, difference in motif score count, 
	indication if the variant falls in a high information position, 
	difference in motif score (between the variant and ref), 
	the overlapping domains (between the variant and ref), reference AA, chancged AA, 
	grantham score, SIFTcat, SIFTval, cDNA position, , CDS position, protein position. 

Abbreviation used for the consequences are:
	SG, Stop_Gained;				NS, Non_Synonymous (missense);
	IF, Inframe_Insertion;			FS, Frame_Shift;
	SL, Stop_Lost;					CS, Canonical_Splice (splice donor);
	S, Splice_Site (splice donor);	NC, Noncoding_Change (non-coding exon);
	SN, Synonymous;					IG, Intergenic;
	DN, Downstream;					UP, Upstream;
	R, Regulatory_Region;			U5, 5Prime_UTR;
	U3, 3Prime_UTR;					I, Intronic;
	O, Unknown. 
	

:Edited by: Seyan Hu
:Date: 29-11-2022
:Usage: python <script>.py -v <VEP output> -r <Reference chromosome> -g <grantham file>

"""
# Import dependencies. 
import sys, os
from optparse import OptionParser
from collections import defaultdict
import pysam
import string
import random
from itertools import tee
import pdb
import numpy as np
import copy


# OptionParser for input files. 
parser = OptionParser("%prog [options]")
parser.add_option("-v", "--vep", dest="vep", help="File with variant effect predictor output (def '')",default="") 
parser.add_option("-r","--reference", dest="reference", help="Path to reference sequences",default='')
parser.add_option("-g", "--grantham", dest="grantham", help="Path to Grantham score annotation file",default='')

parser.add_option("-a","--all", dest="all", help="Produce all output lines, rather than random one on same hierarchy (def OFF)",default=False,action="store_true")
parser.add_option("-b","--verbose", dest="verbose", help="Turn verbose messages on (def OFF)",default=False,action="store_true")

(options, args) = parser.parse_args()

print("Starting on VEP processing", flush=True)

# Defines path to VEP output. 
path_vep = options.vep
print("ref: " + options.reference, flush=True)

# Define chromosome. 
chrom = options.reference.split('ref')[-1]
chrom = chrom.split('.')[0]
print(chrom, flush=True)

# Define headers for output file. (Originally there were also 'Dst2Splice', 'Dst2SplType', these were removed since this information is already in the consequences of the VEP output)
elist = ['#Chrom','Pos', 'Ref', 'Alt', 'isTv', 'Consequence', 'GC', 'CpG', 'motifECount', 'motifEHIPos', 'motifEScoreChng', 'Domain', 'oAA','nAA', 'Grantham','SIFTcat', 'SIFTval','cDNApos', 'relcDNApos', 'CDSpos','relCDSpos', 'protPos', 'relprotPos']


# List for transversions and transitions. 
transversions = set([('A','C'),('C','A'),('T','A'),('A','T'),('C','G'),('G','C'),('G','T'),('T','G')])
transitions = set([('C','T'),('T','C'),('G','A'),('A','G')])


# List of hierachy of the consequences. 
hierachy1 = ["SG","CS","NS","SN","FS","SL","S","IF","U5","U3","R","IG","NC","I","UP","DN","O"]


# Function for counting GC and CpG sites in a window of 75 bases (75 bases before and after the variant).
# Returns percentage GC and CpG counts inside this window. 
def count_GC_CpG(chrom,start,end,previous_pos,previous_GC,previous_CpG,window,seq_tabix):
	if previous_pos != start:
		try:
			sequence = seq_tabix.fetch(chrom,start-window,end+window)
			CpG,GC = 0,0
			count = 0
			lbase = ''
			for base in sequence:
				count += 1
				if base in 'GC': GC+=1
				if lbase == 'C' and base == 'G': CpG+=1
				lbase = base
			
			if count > 0: 
				return GC/float(count),CpG/(count*0.5)
			else: 
				return '-','-'
		except:
			return '-','-'
	else:
		return previous_GC,previous_CpG


# Function for formatting the header or annotations of a variant and returns a formatted line. 
def annotation2line(data,header=False):
	global elist
	if header:
		return "\t".join(elist)+"\n"
	else:
		fstr = ""
		for elem in elist:
			if elem in data: 
				fstr += "%s\t"%(data[elem])
			else: 
				fstr += "-\t"
	return fstr[:-1] + "\n"


# Function for extracting chromosome, chromosome position and the Ref and Alt alleles from the vcf file and VEP output. 
# Appends data to the given dict. 
# Returns the given dict containing these data. 
def extract_alleles_locs(output_dict,fVCoord,fVallele,chrom,previous_pos,previous_ref_allele,fVName,vepfields): 
	output_dict['#Chrom'] = chrom
	output_dict['Pos'] = int(vepfields[fVCoord].split(':')[1])
	output_dict['Alt'] = vepfields[fVallele].upper()
	
	#if position has changed, change reference, otherwise keep the previous one.
	if previous_pos != output_dict['Pos']:
		output_dict['Ref'] = vepfields[fVName].split('/')[0][-1]
	else:
		output_dict['Ref'] = previous_ref_allele


# Function for extracting features from VEP annotated file with the given labels. 
# Appends features to given dict.
# Returns the dict with the appended features. 
def extract_transcript_coding_prot_feature(output_dict,vepfields,position,label1,label2,pos,previous_pos,previous_label1,previous_label2):
	if vepfields[position][0]!="-":
		helper = []
		elength = None
		helper = [x.strip() for x in vepfields[position].split("/")]
		if len(helper) == 2:
			elength = int(helper[-1])
			vepfields[position] = helper[0]
		else:
			sys.exit('At this point a list called annotations should be iterated but this is not defined in this script or in the original. Chrom:%s Pos:%s :%s cDNA%s' %(output_dict['#Chrom'],output_dict['Pos'],vepfields[position]))
		output_dict[label1] = vepfields[position].replace('?-','').replace('-?','').split('-')[0]
		if elength != None:
			output_dict[label2] = "%.2f"%(min(1.0,float(vepfields[position].replace('?-','').replace('-?','').split('-')[0])/elength))
		
		return output_dict
	else:
		output_dict[label1] = "-"
		output_dict[label2] = "-"
		return output_dict


# Function for Extracting the consequences from the VEP annotated vcf file.
# It also gives the consequence an abbreviation. 
# Appends to the given dict these consequences and returns it. 
def extract_consequences(output_dict,vepfields,fVconseq):
	consequences = set([x.strip() for x in vepfields[fVconseq].split(",")])
	if ("stop_gained" in consequences) or ("start_lost" in consequences):
		output_dict["Consequence"] = "SG"
	elif ("missense_variant" in consequences) or ("initiator_codon_variant" in consequences) or ("protein_altering_variant" in consequences):
		output_dict["Consequence"] = "NS"
	elif ("inframe_insertion" in consequences) or ("inframe_deletion" in consequences):
		output_dict["Consequence"] = "IF"
	elif ("frameshift_variant" in consequences):
		output_dict["Consequence"] = "FS"
	elif ("stop_lost" in consequences) or "incomplete_terminal_codon_variant" in consequences:
		output_dict["Consequence"] = "SL"
	elif ("splice_donor_variant" in consequences) or ("splice_acceptor_variant" in consequences):
		output_dict["Consequence"] = "CS"
	elif ("splice_region_variant" in consequences):
		output_dict["Consequence"] = "S"
	elif ("non_coding_exon_variant" in consequences) or ("mature_miRNA_variant" in consequences) or ("non_coding_transcript_exon_variant" in consequences):
		output_dict["Consequence"] = "NC"
	elif ("synonymous_variant" in consequences) or ("stop_retained_variant" in consequences):
		output_dict["Consequence"] = "SN"
	elif ('intergenic_variant' in consequences ):
		output_dict["Consequence"] = "IG"
	elif ('downstream_gene_variant' in consequences):
		output_dict["Consequence"] = "DN"
	elif ('upstream_gene_variant' in consequences):
		output_dict["Consequence"] = "UP"
	elif (("feature_truncation" in consequences) or ("feature_elongation" in consequences)):
		if ("coding_sequence_variant" in consequences):
			output_dict["Consequence"] = "O"
		elif ("regulatory_region_variant" in consequences) or ("TF_binding_site_variant" in consequences) or ("regulatory_region_amplification" in consequences) or ("feature_elongation" in consequences):
			output_dict["Consequence"] = "R"
		elif "5_prime_UTR_variant" in consequences:
			output_dict["Consequence"] = "U5"
		elif "3_prime_UTR_variant" in consequences:
			output_dict["Consequence"] = "U3"
		elif "intron_variant" in consequences:
			output_dict["Consequence"] = "I"
		elif ("non_coding_transcript_variant" in consequences):
			output_dict["Consequence"] = "NC"
		elif ('intergenic_variant' in consequences ):
			output_dict["Consequence"] = "IG"
		elif ('downstream_gene_variant' in consequences):
			output_dict["Consequence"] = "DN"
		elif ('upstream_gene_variant' in consequences):
			output_dict["Consequence"] = "UP"
		else:
			output_dict["Consequence"] = [x.strip() for x in vepline[fVconseq].split('_')[0].upper()]
			sys.exit('This variant has an unrecognized Consequence, %s '%('\t'.join(output_dict)))
	elif ("regulatory_region_variant" in consequences) or ("TF_binding_site_variant" in consequences) or ("regulatory_region_amplification" in consequences) or ("feature_elongation" in consequences):
		output_dict["Consequence"] = "R"
	elif "5_prime_UTR_variant" in consequences:
		output_dict["Consequence"] = "U5"
	elif "3_prime_UTR_variant" in consequences:
		output_dict["Consequence"] = "U3"
	elif "intron_variant" in consequences:
		output_dict["Consequence"] = "I"
	elif ("coding_sequence_variant" in consequences):
		output_dict["Consequence"] = "O"
	elif ("non_coding_transcript_variant" in consequences):
		output_dict["Consequence"] = "NC"
	elif ('intergenic_variant' in consequences ):
		output_dict["Consequence"] = "IG"
	elif ('downstream_gene_variant' in consequences):
		output_dict["Consequence"] = "DN"
	elif ('upstream_gene_variant' in consequences):
		output_dict["Consequence"] = "UP"
		
	elif len(consequences) == 1:
		output_dict["Consequence"] = [x.strip() for x in vepline[fVconseq].split('_')[0].upper()]
		sys.exit('This variant has an unrecognized Consequence, '%(output_dict))
	else:
		sys.stderr.write("Need simplification: %s %s)\n"%(consequences,vepline[fVconseq]))
		output_dict["Consequence"] = "O"
	return output_dict


# Function for extracting the changed amino acid from VEP. 
# Appends the AA before and after the mutation to the given dict and returns it. 
def extract_Aminoacids(output_dict,vepfields,fVAA,previous_pos,previous_ref,previous_alt,previous_oAA,previous_nAA):
	if vepfields[fVAA] != "-": 
		hfields = [x.strip() for x in vepfields[fVAA].split('/')]
		if len(hfields) == 2:
			output_dict['oAA'] = hfields[0]
			output_dict['nAA'] = hfields[1]
		elif hfields[0] != "X":
			output_dict['oAA'] = hfields[0]
			output_dict['nAA'] = hfields[0]
		elif len(hfields[0]) == 1:
			output_dict['oAA'],output_dict['nAA'] = ("-","-")
		else:
			print(hfields,output_dict)
			sys.exit('The field resevered for Aminoacid has a corrupted value. The value and the latest output dict are printed.')
	else:
		output_dict['oAA'],output_dict['nAA'] = ("-","-")
	
	return output_dict


# Function for extracting the extra information in the VEP annotation. 
# It appends the extra info to the given dict and returns it. 
def extract_extra(output_dict,vepfields,fVExtra):
	
	for elem in [x.strip() for x in vepfields[fVExtra].split(";")]:
		hfields = [x.strip() for x in elem.split("=")]
		
		if hfields[0] == "SIFT":
			hfields2 = [x.strip() for x in hfields[1].rstrip(")").split("(")]
			output_dict['SIFTcat'] = hfields2[0]
			output_dict['SIFTval'] = hfields2[1]
		elif hfields[0] == "DOMAINS":
			ncoils,tmhmm,sigp,ndomain,lcompl = False,False,False,False,False
			
			for dfields in [x.strip() for x in hfields[1].split(',')]:
				if len([x.strip() for x in dfields.split(":")]) < 2:
					continue
				
				category,name = [x.strip() for x in dfields.split(":")][0:2]
				if ("_domain" in category) or ("_profile" in category):
					ndomain = True
				else:
					name = name.lower()
					if "coil" in name: ncoils = True
					elif "tmhelix" in name: tmhmm = True
					elif "signalp" in name: sigp = True
					elif name == "seg": lcompl = True
			
			# Implement simple hierarchy of domain annotations:
			if ncoils: output_dict['Domain'] = "ncoils"
			elif tmhmm: output_dict['Domain'] = "tmhmm"
			elif sigp: output_dict['Domain'] = "sigp"
			elif ndomain: output_dict['Domain'] = "ndomain"
			elif lcompl: output_dict['Domain'] = "lcompl"
		
		elif hfields[0] == "HIGH_INF_POS":
			output_dict['motifEHIPos'] = "True" if "Y" == hfields[1] else "False"
		elif hfields[0] == "MOTIF_SCORE_CHANGE":
			output_dict['motifEScoreChng'] = hfields[1]
			output_dict['motifECount'] = '1'
	
	#if Extras are not defined, "-" has to be stored
	if 'motifEScoreChng' not in output_dict:
		output_dict['motifEScoreChng'] = "0.0"
		output_dict['motifECount'] = '0'
	if 'motifEHIPos' not in output_dict:
		output_dict['motifEHIPos'] = "-"
	
	if 'Domain' not in output_dict:
		output_dict['Domain'] = "-"
	if 'SIFTval' not in output_dict:
		output_dict['SIFTval'] = "-"
	if 'SIFTcat' not in output_dict:
		output_dict['SIFTcat'] = "-"
	return output_dict
	

# Function for returning the most deleterious annotation for the same variant, 
# when there are two annotations given for a single variant. 
def indexing(previous,current):
	global hierachy1
	
	index_current = hierachy1.index(current)
	index_previous = hierachy1.index(previous) 
	
	if index_previous > index_current:
		return current
	else:
		return previous


# Function for Reading the grantham file and returns it as a dict. 
def read_grantham(filename):
	grantham = {}
	if os.path.exists(filename):
		infile = open(filename)
		for line in infile:
			fields = line.split()
			if len(fields) == 2:
				AAs = tuple(fields[0].upper().split("-"))
				grantham[AAs]=fields[1]
			else:
				sys.stderr.write("Grantham scores, unexpected line, skipping: %s\n"%line.strip())
		infile.close()
	else:
		sys.stderr.write("Grantham scores input file does not exist: %s\n"%filename)
	return grantham


# Open reference of chr, generated variants (vcf without annotations), exome files. 
ref_fasta = pysam.Fastafile(options.reference)
print("Read reference fasta", flush=True)

# Reads in the grantham matrix.
grantham_matrix = read_grantham(options.grantham)
print("Read grantham", flush=True)

# Open VEP annotated variants. 
# Header VEP: #Chrom	Start	End	Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
vepinput = open(path_vep, 'r')
fline = True
print("Read VEP annotated variants", flush=True)

#processes the body of the VEP file and writes out the annotations.
previous_dict = {'#Chrom':'-', 'Pos':0, 'Ref':'X', 'Alt':'X', 'isTv':'-', 'Consequence':'-', 'GC':'-', 'CpG':'-', 'motifECount':'-', 'motifEHIPos':'-', 'motifEScoreChng':'-', 'Domain':'-', 'oAA':'-', 'nAA':'-', 'Grantham':'-', 'SIFTcat':'-', 'SIFTval':'-', 'cDNApos':'-', 'relcDNApos':'-', 'CDSpos':'-', 'relCDSpos':'-', 'protPos':'-','relprotPos':'-' }
previous_minscore = None


#I have found a variant which causes problems because it has more annotated fields than the others. I will sort them out and treat them afterwards
problem_out = open('problem_out.tsv','a')


# Create output. 
outp_chr = options.reference.split('ref')[-1]
outp_chr = outp_chr.split('.')[0]
print(outp_chr, flush=True)
outp = open(outp_chr + '_vep_processed.vcf', 'w')


#processes the VEP header
for vepline in vepinput:
	if vepline.startswith('##'): 
		continue
	
	# Was previously 'Uploaded_variation' instead of Chrom. 
	if vepline.startswith('#Chrom') & fline: #treats prints first line 
		vepfields = [x.strip() for x in vepline.split('\t')]
		
		# Uploaded_var location in the fields set to 3, was previously 0 (+3 to everything). 
		# Length of field in vep was previously set to 14. 
		if len(vepfields) == 17:
			
			fVName = 3				# Uploaded_variation 
			fVCoord = 4	 			# Location variation
			fVallele = 5			# Allele
			fVgene = 6				# Gene
			fVfeature = 7 			# Feature
			fVfeattype= 8 			# Feature_type
			fVconseq = 9			# Consequence
			fVcDNA = 10	 			# cDNA_position
			fVCDSPos = 11 			# CDS_position
			fVpPOS = 12	 			# Protein_position 
			fVAA = 13				# Amino_acids
			fVCodon = 14			# Codons
			fVVar = 15				# Existing_variation
			fVExtra = 16			# Extra
		else:
			sys.exit("Unknown input data format.")
		fline = False
		
		# Writes header to standard output.
		outp.write(annotation2line({},True)) 
		break



fline = True
# Iterates over the lines in the VEP output. 
for vepline in vepinput:
	
	# Splits the line into fields. 
	vepfields = [x.strip() for x in vepline.split('\t')]

	
	# Check if vep output has correct number of fields. 
	if len(vepfields) != 17:
		problem_out.write(vepline)
		continue
		
		
	# Here begins the parsing of the vep input line, the function consists of many subfunctions for each feature respectively. 
	output_dict = {}
	# Every time when previous position and current one are different, go one further on the reference allele
	
	
	# Performs the function for extracting chromosome, chromosome position and the Ref and Alt alleles from the vcf file and VEP output.
	extract_alleles_locs(output_dict,fVCoord,fVallele,chrom,previous_dict['Pos'],previous_dict['Ref'],fVName,vepfields)
	
	
	# Performs function that returns the percentage of GC and CpG in given window of 75 (75 bases before and after the position of the variant).
	# Replaced 'ref_ident' with output_dict['#Chrom'], otherwise it would not work if the reference is from ncbi.
	output_dict['GC'],output_dict['CpG'] = count_GC_CpG(output_dict['#Chrom'],output_dict['Pos'],output_dict['Pos'],previous_dict['Pos'],previous_dict['GC'],previous_dict['CpG'],75,ref_fasta)
	
	
	# Checks whenever the mutations is a transversion or not. 
	output_dict['isTv'] = (output_dict['Ref'],output_dict['Alt']) in transversions
	
	
	# Performs function that returns features from VEP annotated file with the given labels.
	output_dict = extract_transcript_coding_prot_feature(output_dict,vepfields,fVcDNA,'cDNApos','relcDNApos',output_dict['Pos'],previous_dict['Pos'],previous_dict['cDNApos'],previous_dict['relcDNApos'])
	output_dict = extract_transcript_coding_prot_feature(output_dict,vepfields,fVCDSPos,'CDSpos','relCDSpos',output_dict['Pos'],previous_dict['Pos'],previous_dict['CDSpos'],previous_dict['relCDSpos'])
	output_dict = extract_transcript_coding_prot_feature(output_dict,vepfields,fVpPOS,'protPos', 'relprotPos',output_dict['Pos'],previous_dict['Pos'],previous_dict['protPos'],previous_dict['relprotPos'])
	
	
	# Performs function that returns the AA before and after the mutation.
	output_dict = extract_Aminoacids(output_dict,vepfields,fVAA,previous_dict['Pos'],previous_dict['Ref'],previous_dict['Alt'],previous_dict['oAA'],previous_dict['nAA'])
	
	
	# Appends the correct Grantham score to the change in AA.
	if (output_dict['nAA'],output_dict['oAA']) in grantham_matrix:
		output_dict['Grantham'] = grantham_matrix[output_dict['nAA'],output_dict['oAA']]
	else:
		output_dict['Grantham'] = '-'
	
	
	# Performs the function for extracting the consequences from the VEP output.
	output_dict = extract_consequences(output_dict,vepfields,fVconseq)
	
	# Performs the function to process the extra field in VEP output. 
	output_dict = extract_extra(output_dict,vepfields,fVExtra)
	
	
	# If it is the first line to be parsed, store it in previous_dict and move on.
	if fline:
		previous_dict = copy.deepcopy(output_dict)
		fline = False
		continue
		
		
	# Checks for duplicates (variants with multiple annotations) 
	# and stores the adjusted line ('output_dict') to 'previous_dict'. 
	elif (previous_dict['Pos'] == output_dict['Pos']) and (previous_dict['Alt'] == output_dict['Alt']):
		
		
		# If statements used regarding other features that may be influenced by duplications such as cDNApos, CDSpos, protPos, oAA, nAA, grantham.
		if (float(previous_dict['motifEScoreChng'])) < (float(output_dict['motifEScoreChng'])):
			output_dict['motifEScoreChng'] = previous_dict['motifEScoreChng']
		
		
		# If duplicate, add up all motifs for the final motifEScore.
		output_dict['motifECount'] = str(int(previous_dict['motifECount']) + int(output_dict['motifECount']))
		
		
		# Performs function that takes the most deleterious VEP Consequence if there are variants with multiple annotations. 
		output_dict['Consequence'] = indexing(previous_dict['Consequence'],output_dict['Consequence'])
		
		for key in previous_dict.keys():
			if output_dict[key] == '-' and previous_dict[key] != '-':
				output_dict[key] = previous_dict[key]
		
		previous_dict = copy.deepcopy(output_dict)
		continue
		
	# If no duplicate output previous_dict and store output_dict in previous_dict
	else: 
		outp.write(annotation2line(previous_dict))
		previous_dict = copy.deepcopy(output_dict)


# Writes last line which is now stored in previous_dict to output. 
outp.write(annotation2line(previous_dict))
vepinput.close()
