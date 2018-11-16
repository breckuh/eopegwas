# preeclampsiagwas
Maternal Cardiovascular-Related Single Nucleotide Polymorphisms, Genes and Pathways Associated with Early-Onset Preeclampsia

## Results

You can download the results of the scripts in the `src` folder here: https://raw.githubusercontent.com/breckuh/preeclampsiagwas/master/src/main.nb.html

Results of snpQC report are here: https://github.com/breckuh/preeclampsiagwas/blob/master/snpQCreport.pdf

## Data required

You'll need these before the code can be run:

- data/072312_FinalReport.txt
- data/clinical.csv
- data/cardio-metabo_chip_11395247_a.csv

## Create some folders:
- ignore/
- cache/
- output/

## Pvals
For some of the R notebooks, you'll first have to run `thread.r` to generate the pvals for each model you want to test. `thread.r` will save
results to the cache folder.

## Phylotree
You'll need to download and install this program: http://chibba.pgml.uga.edu/snphylo/ to run the phylotree.rmd script. Put the Phylotree scripts in the ignore folder.

## Testing
Run `rmarkdown::render('main.rmd')` to execute the majority of the code