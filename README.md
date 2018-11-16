# preeclampsiagwas
Maternal Cardiovascular-Related Single Nucleotide Polymorphisms, Genes and Pathways Associated with Early-Onset Preeclampsia

## Data required

You'll need these before the code can be run:

- data/072312_FinalReport.txt
- data/clinical.csv
- 

## Create some folders:
- ignore/
- cache/
- output/

## Pvals
For some notebooks, you'll also first have to run `thread.r` to generate the pvals for each model you want to test. That script will save
results to the cache folder.

## Phylotree
You'll need to download and install this program: http://chibba.pgml.uga.edu/snphylo/ to run the phylotree.rmd script. Put it in the ignore folder.

## Testing
Run `rmarkdown::render('main.rmd')` to execute the majority of the code