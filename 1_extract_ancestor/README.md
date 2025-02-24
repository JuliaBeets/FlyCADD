# Extract reconstructed ancestral sequence from multi-species sequence alignment

The first step in creating the FlyCADD was obtaining the reconstructed ancestral sequence. The scripts in this folder were used to extract the ancestral sequence from the 166-way Multi-species Sequence Alignment (MSA) provided as an early version of the 298-way MSA by B. Kim et al. (2024). The alignment was generated with Cactus and includes 166 Drosophila species and a reconstructed ancestral sequence for each node in the phylogeny, all in `.maf` format. Species names used in these scripts should be written as they are in the MSA.

### Mark outgroup
The MSA that was used contained reconstructed ancestral sequences. The ancestral node of interest was marked with `mark_outpgroup.py` to make this sequence recognizable for the next steps. The input for the script is the MSA in `.maf` format and the output is the MAF-file with one ancestor marked by the prefix "Ancestor_".  <br />
Usage: `python mark_outgroup.py -a <name of ancestor in MSA> -f <reference species name in MSA> -p <path to files> -i <.maf file, e.g. drosophila.maf.gz>`

### Apply mafTools
The marked `.maf` file was preprocessed with three tools from dentearl/mafTools (doi: 10.1101/gr.174920.114):
- mafDuplicateFilter: filter out duplications
- mafStrander: ensure alignment blocks are relative to the positive strand of the reference species
- mafRowOrderer: order alignment blocks within the MSA according to the order provided

The input of the script was the marked MSA in `.maf` format from the previous step and the output is a processed `.maf` file with filtered and ordered alignment blocks. <br />
Usage: `python apply_mafTools.py -m <path to marked MSA> -g <reference species name> -o <order of species for RowOrderer: ref_species, ancestor, e.g. D_MELANOGASTER,Ancestor> -s <path to file after mafStrander, e.g. ./${anc}/mSTR/mSTR_> -r <path to file after mafRowOrderer, e.g. ./${anc}/mRO/mRO_> -c <clean yes or no> -f <path to file after mafDuplicateFilter, e.g. ./${anc}/mDF/mDF_>`

### Sort alignment blocks by chromosome
The preprocessed `.maf` file from the previous step was sorted based on the chromosome order of the reference species of interest. The input for the script is the preprocessed `.maf` file from the previous step and the output are separate `.maf` files per reference species chromosome.<br />
Usage: `python sort_by_chromosome.py -p <path to preprocessed MSA from the previous step> -f <prefix of preprocessed MSA from the previous step> -s <name of reference species of interest> -c <list of sorted chromosomes, e.g. '2L,2R,3L,3R,4,X,Y'> -o <path for output> -c <clean yes or no>`

### Sort alignment blocks by genomic coordinates
The contents of the `.maf` files from the previous step were reordered based on the genomic positions of the reference species of interest. The input are `.maf` files per chromosome and the output files are `.maf` files per chromosome with ordered contents based on genomic positions.<br />
Usage: `python sort_msa_blocks.py -p <path to preprocessed MSA from the previous step> -s <name of reference species> -o <path for output>`

### Remove unwanted species
This script removed all unwanted species from the `.maf` files to make sure only the reference species and ancestor sequences were left in the `.maf` files.<br />
Usage:` python remove_species.py -p <path to preprocessed MSA from the previous step> -f <prefix of preprocessed MSA from the previous step> -s <name of reference species> -r <path for output>`

### Remove reverse strand of reference species
This script removed alignment blocks that were based on the reverse strand of the reference species. The output of the script is the final processed MSA in `.maf` format with the alignment blocks including the forward strand of the reference species and the aligned ancestral sequence of the ancestor of interest. These files were split per chromosome, sorted based on genomic location and filtered.<br />
Usage: `python remove_opposite_strand.py -p <path to preprocessed MSA from the previous step> -f <file prefix of preprocessed MSA from the previous step> -r <path for output>`

### Extract ancestral sequence
The `extract_ancestor.py` script extracted the ancestral sequences based on the position of the reference species they mapped to and resulted in FASTA files per chromosome with the reconstructed ancestral sequence. Genomic positions are relative to the coordinates of the reference species (dm6). To extract the reconstructed ancestor for each chromosome, a wrapper for this script ensures the script is executed for each chromosome.<br />
Usage: `python wrapper_extract_ancestor.py -p <path to processed MSA files from the previous step> -a <name of ancestor, e.g. Ancestor_Anc117> -s <name of reference species> -f <prefix of processed MSA files from the previous step> -g <path to wrapped python script> -i <path to the FASTA-index of the reference genome>`