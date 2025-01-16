#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 14-10-2022
:Usage: python <script.py> <vcf.gz file> <chr list>

Creates for each chr a frequency file from a vcf file.

:Example:
python scripts/gen_freq.py data/vcf_data/mgp_REL2021_snps.vcf.gz '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y'
'''
# Import dependencies
import sys, os
from optparse import OptionParser

# OptionParser for the directories of the input 
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="Path to individual level vcf file", default = "reference.vcf")
parser.add_option("-c", "--chromosomes", dest="chromosomes", help="List of chromosomes to be considered when extracting freq data", default = '1,2,3,4,5')

(options, args) = parser.parse_args()


# Creates a list for the chromosomes
chr_list = options.chromosomes.split(',')


# Loop through list of chr and perform vcftools
for chr_num in chr_list:
	os.system('vcftools --gzvcf ' + options.vcf + ' --chr ' + chr_num + ' --remove-indels --non-ref-af 0.9 --max-non-ref-af 1.0 --stdout --freq > freq_' + chr_num + '.out')
	os.system('vcftools --gzvcf ' + options.vcf + ' --chr ' + chr_num + ' --remove-indels --non-ref-af 0.0 --max-non-ref-af 0.1 --stdout --freq > freq_reversed_' + chr_num + '.out')