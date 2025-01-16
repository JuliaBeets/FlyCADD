# Generate genome-wide CADD-scores for all possible SNVs
The fifth and final step in creating FlyCADD is applying the trained ML model to all possible variants of the reference genome and computing final CADD-scores. All variants of the reference genome are generated and annotated with the same annotations as the derived and simulated variants (the training and testing datasets) and the ML model is used to predict the probability for the variant to be in category 0 (deleterious mutations).

### Generate all possible variants of the reference genome
For each position on the reference genome, all three possible variants are generated on chromosome 2L, 2R, 3L, 3R, 4 and X. The script generates all variants based on the reference genome (per chromosome). <br/>
Usage: `python gen_all_variants.py -i <path to directory with reference genome>`

### VEP annotation
Each variant is annotated with VEP and with the second script the output is processed to return a tab-delimited, encoded file with additional annotations such as Grantham score for the effect of substitutions between amino acids based on chemical properties (obtained from DOI: 10.1126/science.185.4154.862). The final list of annotations at this stage is: chromosome, position, reference allele, alternative allele, is it transversion, GC%, CpG%, difference in motif score (bound by TF), is it high information position, overlapping domains (between the variant and reference), reference AA, variant AA, grantham score, cDNA position, , CDS position, protein position. <br/>
Usage: `python vep_wrapper.py -i <path to directory with all generated variants (VCF)> -s <species name, e.g. drosophila_melanogaster>`
Usage: `python vep_output_processing_wrapper.py -v <path to directory with VEP-annotated VCF files> -r <path to directory with reference genome> -g <path to Grantham table .tsv format>`

### Merge all annotations
For each variant in the VEP annotated, processed VCF files the PhastCons, PhyloP, repeat positions and gene positions (generated at step 3) are added. This script can be expanded to add more annotations at a later stage. <br/>
Usage: `python merge_all_annotations.py -v <path to directory with VEP processed VCF files> -p <path to directory with processed PhastCons and PhyloP scores> -r <path to directory with repeat positions> -g <path to directory with gene positions>`

### Add annotations
Adds repeattype, TFBS, miRNA, ReMap, CRM, 124Conservation, BG3 and S2 chromatin states and PhyloP and GERP computed from MSA.
Usage: python add_annotations.py -o <Path to merged annotations> -a <Path to annotation to add> -b <Path to second annotation to add> -c <Path to third annotation to add> -d <Path to fourth annotation to add> etc...

### Split dataset in equal chunks
To be able to run next steps in parallel, the VCF files are split into chunks of equal size (1.000.000 variants). <br/>
Usage: `python chunk_DF.py -i <path to directory containing fully annotated VCF files> -l <number of lines per chunk>`

### Scaling and imputation of all variant annotations
Based on the mean of the simulated variants, missing values are handled with imputation and categorical features with one-hot encoding. This script is wrapped to process all files.<br/> 
Usage: `python encoding_mv_wrapper.py -i <path to chunked dataframes> & ... ` continue for as many directories with chunks exist (depending on the number of lines per chunk)<br/>
The features are scaled with the stdev used previously for the scaling of the derived and simulated data. Empty features (that only contains 0) are removed from the DF. <br/>
Usage: `python feature_scaling.py -i <path to first directory with chunks> -d <path to directory with model estimates after training> & ...` continue for as many directories with chunks exist (depending on the number of lines per chunk)

### Prediction of posterior likelihood for being a deleterious variant by the trained ML model
The previously trained model is loaded and applied to the fully annotated and processed whole genome variants to predict the posterior likelihood of the variants being neutral or benign (derived, class 1) or deleterious (simulated, class 0). The output file consists of the columns probability for class 0, probability for class 1, position, reference allele, alternative allele.
Usage: per folder `python MLpredict.py -i <path to scaled variants> -m <path to model file>

The predictions are being processed to remove the predictions for class 1 and keep the probability of class 0 (deleterious).<br/>
Usage: python process_prediction.py -i <path to initial predictions>

### Merge all variants
This script merges the chunks per chromosome.<br/>
Usage: `python merge_chunks.py -i <Path to processed chunks> 