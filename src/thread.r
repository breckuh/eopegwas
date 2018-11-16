#!/usr/bin/env Rscript
source("setup.r")

# We have this tiny script because SNPassoc takes hours per model when you have >100k SNPs

args = commandArgs(trailingOnly=TRUE)
modelName = args[1]
cacheKey = args[1]

# ./thread.r codominant
# ./thread.r dominant
# ./thread.r recessive
# ./thread.r log-additive
# ./thread.r overdominant


# snpTable = pe.getSnpTableFilip(callRate = .95)
snpTable = pe.getSnpTable(callRate = .95)

cache.get(glue("{cacheKey}Filip"), function () SNPassoc::WGassociation(preeclampsia, data = snpTable, model = modelName))
