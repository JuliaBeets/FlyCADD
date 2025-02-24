# Scripts FlyCADD 

This repository contains the scripts, tools and the workflow used to create the FlyCADD model impact prediction of single nucleotide variants in *Drosophila melanogaster*. This GitHub repository contains the directions for FlyCADD development. A locally excecutable pipeline for scoring novel variants with FlyCADD can be found on Zenodo 10.5281/zenodo.14887337. 

For the whole pipeline, an unmasked reference genome and an alignment file in `.maf` format were used, except when specified otherwise. All annotations and FlyCADD scores are 1-based.
______________________________________________________________________________________________________________________
The repository represents FlyCADD as published in the original manuscript: {will be added once published}. </br>
The precomputed FlyCADD scores, alignment file and other supporting data can be found on Zenodo 10.5281/zenodo.14887337. 

Contact: Julia Beets - j.beets@vu.nl 
______________________________________________________________________________________________________________________

## How to use

FlyCADD development is described in {DOI added once published}, refer to the manuscript for details. The five steps for creating FlyCADD scores include:
- [Extract ancestor](1_extract_ancestor/): The ancestral sequence was extracted from the Multiple Alignment Format file of a 166-way multi-species alignment (Kim et al., 2024). Each node in the alignment file can be chosen, for FlyCADD we extracted the ancestral sequence at node 117, the last common ancestor between *D. melanogaster* and *D. tani*. 
- [Create trainingdata](2_create_trainingdata/): The training dataset consists of derived and simulated variants. Derived variants were identified as those with a different allele in the ancestral sequence, at high frequency (> 0.9) in *D. melanogaster* populations. Simulated variants were simulated *de novo* based on mutation rates extracted by comparison of the reference genome dm6 and reconstructed ancestral genome. The set of simulated variants was randomly trimmed to match the number of derived variants.
- [Annotate variants](3_annotate/): The derived and simulated variants were annotated with 691 features, consisting of 38 distinct annotations and combinations thereof.
- [Train and test model](4_train_test_model/): 90% of the derived and simulated variants were used for training the model. The model is a logistic regressor by Turi Create. The remaining 10% of the dataset was used for testing the model using the built-in function model.evaluate().
- [Generate genome-wide FlyCADD scores](5_generate_genome-wide_FlyCADDscores/): Genome-wide FlyCADD scores were generated by determining all possible variants on the reference genome dm6, annotating these variants with the same annotations as the training dataset and applying the model to provide a probability of each variant belonging to the simulated class, and thus being impactful.

Each folder includes the scripts and descriptions. A [`.yml`](FlyCADD-environment.yml) file provides the user with the Conda-environment containing the packages required for FlyCADD development or excecution.
These scripts were excecuted on an High-Performance Cluster (Linux) in the order and with the parameters described in each respective folder.

## How to cite

{will be added once published} 

## Licence
The usage of scripts of CADD (Kircher et al. 2014) is indicated in the script descriptions. The licence is in LICENCE_hCADD.

Original FlyCADD scripts are licenced in LICENCE_FlyCADD.