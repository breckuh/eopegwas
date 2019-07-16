# preeclampsiagwas

This repository holds the code for the paper: "Maternal Cardiovascular-Related Single Nucleotide Polymorphisms, Genes and Pathways Associated with Early-Onset Preeclampsia".

## Results

The paper contains a summary and explanation of the findings and the `src` folder contains the scripts. You can also see the results of the R code, without running it yourself, in notebook format here: [https://breckuh.github.io/eopegwas/src/main.nb.html](https://breckuh.github.io/eopegwas/src/main.nb.html)

Results of snpQC report are here: [https://github.com/breckuh/preeclampsiagwas/blob/master/snpQCreport.pdf](https://github.com/breckuh/preeclampsiagwas/blob/master/snpQCreport.pdf)

## Data required

You can reuse this code with your own data to perform similar experiments, or use it with the provided simulated data (coming soon!) which is identical to our data except simulated for privacy reasons. Or you can contact us for our raw data to get the exact same results. The data sources required are:

- data/072312_FinalReport.txt
- data/clinical.csv
- data/cardio-metabo_chip_11395247_a.csv

The format for the data is: (todo: add).

## Setup Instructions

To run the R code in this repo:

    git clone https://github.com/breckuh/eopegwas
    cd eopegwas
    mkdir ignore
    mkdir cache
    mkdir output
    rmarkdown::render('main.rmd')

## Pvals
For some of the R notebooks, you'll first have to run `thread.r` to generate the pvals for each model you want to test. `thread.r` will save
results to the cache folder.

## Phylotree
You'll need to download and install this program: [http://chibba.pgml.uga.edu/snphylo/](http://chibba.pgml.uga.edu/snphylo/) to run the phylotree.rmd script. Put the Phylotree scripts in the ignore folder.

## Testing
Run `rmarkdown::render('main.rmd')` to execute the majority of the code
