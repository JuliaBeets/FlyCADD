# Annotate derived and simulated variant sets
The third step in creating the FlyCADD model was annotation of the sets of derived and simulated variants. The variants were annotated with annotations from Ensembl Variant Effect Predictor (VEP), PhastCons, PhyloP, repeats and gene locations, and more. In total, 38 annotations with 691 combined features. All these annotations were processed and merged to create one `.csv` file per chromosome. These annotations with similar scripts are also used for annotation of all possible variants of the reference genome (Step 4).

## VEP annotations
The Variant Effect Predictor by Ensembl annotates variants for which genes and transcripts it affects, the location of the variant (e.g. upstream of a transcript, in coding sequence, in non-coding RNA, in regulatory regions), the consequence of the variation (e.g. stop-gained, synonymous, regulatory region, frame shift, intergenic, no consequence or non-coding change), the allele and more (https://www.ensembl.org/info/docs/tools/vep/index.html). These scripts use VEP to annotate derived and simulated variants per chromosome and process these to create a dataframe where other annotations can be added.
! To use VEP in the offline mode, an offline library for VEP on *D. melanogaster*  was installed from https://ftp.ensembl.org/pub/release-105/variation/vep/drosophila_melanogaster_vep_105_BDGP6.32.tar.gz.

### VEP annotating
This script performs VEP for the derived and simulated variants. The simulated and derived variants (per chromosome) should all be in one directory that is used as input for the script. The output is annotated VCF files containing derived or simulated variants per chromosome.
Usage: `python vep_wrapper.py -i <path to VCF files (derived and simulated) from step 2> -s <species name, e.g. drosophila_melanogaster>`

### VEP output processing
The script `VEP-processing.py` processes the VEP-annotated VCF files per chromosome to return a tab-delimited, encoded file with additional annotations such as Grantham score for the effect of substitutions between amino acids based on chemical properties (obtained from DOI: 10.1126/science.185.4154.862). <br/>
Usage: `python vep_output_processing_wrapper.py -s <path to directory with VCF files> -v <path to directory with VEP annotated VCF files> -r <path to directory with reference genome -g <path to .tsv with Grantham scores> `
<br/>The final list of annotations at this stage is: chromosome, position, reference allele, alternative allele, is it transversion, GC%, CpG%, difference in motif score (bound by TF), is it high information position, overlapping domains (between the variant and reference), reference AA, variant AA, grantham score, cDNA position, CDS position, protein position.

## PhastCons and PhyloP annotations
PhastCons and PhyloP data was downloaded from UCSC Genome Browser database for Multiz Alignment & Conservation (27 Species) at https://hgdownload.soe.ucsc.edu/goldenPath/dm6/ in BigWig format. These scripts reformat the scores to be able to parse them to the processed VEP output.

### Format PhastCons and PhyloP scores
#### BigWig to Wig
BigWigtoWig reformats the BigWig files to Wig files. bigWigToWig is performed twice, for PhastCons and PhyloP seperately. Usage:`bigWigToWig <path to BigWig file> <phastCons|PhyloP>_out.wig`<br/>
#### Reformat to .txt file
These files are further processed to a `.txt` file containing chromosome, position and score. Usage: `python format_pC_pP_wig.py -w <path to converted Wig file>`
#### Split by chromosome
The pre-processed PhastCons and PhyloP files are split to create seperate files per chromosome, performed twice. Usage: `python split_pC_pP_scores.py -f <path to txt file with PhastCons/PhyloP scores> -c <list of chromosomes, e.g. '2L,2R,3L,3R,X,4'>`

## Repeats annotations
Soft-masked reference genome files (GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fa) were obtained from GenBank for *D. melanogaster* and split by chromosome. This script extracts the positions of repeats from these files by extracting the positions of lowercase nucleotides from the FASTA files.<br/>
Usage: `python get_position_repeats.py -r <path to directory containing reference sequence files>`

## Gene location annotations
From FlyBase Annotation Release 6.32 (http://ftp.flybase.net/releases/FB2020_01/dmel_r6.32/fasta/dmel-all-gene-r6.32.fasta.gz), a FASTA file containing sequences of genes was obtained. This FASTA was converted to a `.bed` file to obtain the positions of alleles inside genes. Output was a `.txt` file per chromosome containing a list of positions inside genes (similar format to repeats).

## Merge all annotations
For each variant in the VEP annotated, processed VCF files the PhastCons, PhyloP, repeat positions, DNAshapeR annotations and gene positions are added. This script can be expanded to add more annotations at a later stage.<br/>
Usage: `python merge_all_annotations.py -v <path to directory with VEP processed VCF files> -p <path to directory with processed PhastCons and PhyloP scores> -r <path to directory with repeat positions> -g <path to directory with gene positions>`

## Add annotations
Adds repeattype, TFBS, miRNA, ReMap, CRM, 124Conservation, BG3 and S2 chromatin states and PhyloP and GERP computed from MSA.
Usage: `python add_annotations.py -o <Path to merged annotations> -a <Path to annotation to add> -b <Path to second annotation to add> -c <Path to third annotation to add> -d <Path to fourth annotation to add> etc...`

## Imputation and one-hot encoding
Fully annotated .tsv files for simulated and derived variants are processed with `data_encoding_mv_handling.py` to handle missing values with imputation based on the mean in simulated variants and categorical features are one-hot encoded. This script is wrapped to process all files.
Usage: `python encoding_mv_wrapper.py -i <path to directory containing fully annotated .tsv files`

The final product at this stage consists of `.csv` files per chromosome for simulated and derived variants with merged annotations for all variants, imputed and hot-encoded. Columns depend on the annotations added.