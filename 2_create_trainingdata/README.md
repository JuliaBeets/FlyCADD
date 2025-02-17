# Create derived and simulated variants
The second step in creating the FlyCADD was creating two sets of variants: derived variants (proxy-neutral/benign) and simulated variants (proxy-deleterious). Variants that are derived over evolutionary time between the existence of the ancestor and the modern *D. melanogaster* and present in the modern species at high frequency (> 90 %, fixed/nearly-fixed) were assumed to be neutral or benign. Simulation based on observed mutation rates and other statistics of sequence evolution (obtained by comparing the ancestral sequence to the reference sequence (Release 6)of *D. melanogaster*, included in this repository) results in simulated *de novo* mutations that are assumed to be enriched in deleterious variants as they have not experienced selection depleting deleteriousness.
This part of the pipeline requires the reconstructed ancestral sequences in FASTA format (`.fa`) per chromosome from the previous step, full reference sequence in FASTA format (`.fa`), reference sequence per chromosome in FASTA format (`.fa`) and Python dependencies listed in `FlyCADD-environment.yml`.

## Simulation of *de novo* variants
### Create parameters for simulation
The `create_parameters.py` script computes the number of substitutions and other statistics per chromosome. It is wrapped with a Python wrapper to output a log file per chromosome. <br />
Usage: `python wrapper_create_parameters.py -a <path to reconstructed ancestral sequences> -r < path to reference genome> -c <list of chromosomes, e.g. '2L,2R,3L,3R,4,X'> -p <prefix of reconstructed ancestor files> -r <prefix of reference genome files> -s <path and prefix of log files>`

### Apply parameters to reference genome (simulation)
This script applies the parameters calculated in the previous step to the reference genome to simulate the chosen number of variants (total number of indels + SNPs). The output of this script is a VCF file containing the simulated variants (both indels and SNPs). <br />
Usage: `python apply_parameters.py -p <path to log files> -i <path to reference genome> -o <path for output VCF file> -c <list of chromosomes, e.g. 2L,2R,3L,3R,4,X> -n <number of simulations, e.g. 22000000>`

### Post-processing VCF-file containing simulated variants
Three steps of post-processing of the VCF file containing simulated variants are required to obtain the finalized set of simulated variants with proxy-deleterious SNPs: filter the variants for SNPs and indels, filter these SNPs having a position in the corresponding ancestral sequence and check the substitution rates.<br />
The first script splits the variants and outputs two VCF files: one with SNPs and one with indels.<br />
Usage: `python split_vcf.py -p <path to directory with simulated variants> -i <filename VCF with simulated variants from the previous step>`<br />
The second script filters these SNPs to have a position in the corresponding reconstructed ancestor and removes those that correspond to a gap in the ancestral sequence. <br />
Usage: `python filter_vcf.py -i <path to VCF with simulated indels> -s <path to VCF with simulated SNPs> -a <path to directory containing reconstructed ancestral sequences>` <br />
The third script checks the substitution rates to make sure they follow the same distribution as the substitution rates that were input for the simulation. <br />
Usage: `python check_substitution_rates.py -i <path to complete set of simulated variants> -l <path to log files> -f <path to VCF with filtered, simulated SNPs>`

## Obtain derived variants
There are several scenarios based on which variants are selected as "derived variants":<br />
Scenario 1: The reference and ancestral allele are the same, but the nearly-fixed (> 90 %) variant in the population is different. <br />
Scenario 2: The reference and ancestral allele are different with the reference allele nearly fixed (> 90 %) in the population. <br />
Scenario 3: The reference and ancestral allele are different and no population frequency information for this position. <br />

### Generate allele frequencies from VCF
This script uses a VCF file based on natural populations with variants that occur in these natural populations and extracts allele frequencies for each variant. The output is a frequency file per chromosome containing the location, alleles and their frequency of all high-frequency (> 90 %) variants in the selected natural population. <br />
Usage: `python generate_frequencies.py -v <path to gzipped VCF file> -c <chromosomes of interest, e.g. '2L','2R','3R','3L','4','X'>`

### Select derived variants
The script `derived_var_gen.py` selects based on the criteria derived variants and generates a VCF file as output. The script is wrapped to combine the derived variants of all chromosomes in one VCF file. <br />
Usage: `python wrapper_derived_gen.py <list of chromosomes, e.g. '2L','2R','3L','3R','X','4'> <path to ancestral sequence> <path to reference genome> <path to frequency files>`

## Finalize variants: equal sets of derived and simulated SNPs
Simulated variants are downsampled to match the number of derived variants. Simulated variants are shuffled and the first x variants are selected for the final set of simulated variants to match the number of derived variants.
Usage: `python trim.py -s <path to directory with simulated variants> -d <path to directory with derived variants>` <br />
Finally, a shell-script `proces_sim.sh` was used to sort the trimmed VCF with simulated variants and split the file into VCF files per chromosome.

The final products of this step are two equally sized VCF files with derived single-nucleotide variants and simulated single-nucleotide variants, split per chromosome.


