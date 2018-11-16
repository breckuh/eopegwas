installIfNeeded = function (packages, installFn = install.packages) {
  newPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(newPackages)) installFn(newPackages)
}

installIfNeeded(c("reshape2", "descr", "data.table", "tidyverse", "lattice", "conflicted", "glue", "RColorBrewer", "ggfortify",
  "corrplot", "summarytools", "scales", "skimr", "SNPassoc", "qqman", "BiocManager"))


library(BiocManager)
installIfNeeded(c("GWASTools", "treeio", "ggtree"), BiocManager::install)

library(conflicted)
library(reshape2)
library(descr)
library(data.table)
library(tidyverse)
library(lattice)
library(glue)
library(RColorBrewer)
library(ggfortify)
library(corrplot)
library(summarytools)
library(scales)
library(skimr)
library(qqman)
library(SNPassoc)

library(treeio)
library(ggtree)
library(GWASTools)

source("cache.r")
cache.init(path = "../cache/")

conflict_prefer("filter", "dplyr")

pe.illuminaDataPath = "../data/072312_FinalReport.txt"
pe.illuminaPath = "../data/cardio-metabo_chip_11395247_a.csv"
pe.clinicalDataPath = "../data/clinical.csv"

pe.getSnpData = function (callRate = 0) {
  data = cache.get("snpdata", function () {
    data = data.table::fread(pe.illuminaDataPath, skip=9, header=T, sep="\t", check.names=F, stringsAsFactors=F, data.table=F)
    data$SNP = NULL
    colnames(data) = c("illuminaSnpLocation", "sampleId", "allele1", "allele2", "avgGenCallScore")
    data
    })
  if (callRate > 0)
  return(data[!is.na(data$avgGenCallScore) & data$avgGenCallScore > callRate, ])
  data
}

pe.getSamplesAsRows = function (...) {
  cache.get("samplesAsRows", function () {
    samplesAsRows = pe.getSnpData(...) %>%
      unite(alleles, c(allele1, allele2), sep = "") %>%
      select(-one_of("avgGenCallScore")) %>%
      spread(illuminaSnpLocation, alleles)
  })
}

pe.dropColumnsThatContainNA = function (df) df[ , colSums(is.na(df)) == 0]
pe.dropColumnsThatContain = function (df, searchString) df[, is.na(mapply('charmatch', searchString, df))]
pe.noncalls.dropEntireColumn = function (df) pe.dropColumnsThatContain(df, "-")
pe.dropMonomorphicColumns = function (df) df[,apply(df, 2, function (col) length(unique(col)) != 1)]

pe.getNumericAlleleMatrix = function (...) {
  cache.get("cleanedRows", function () {
    rowsAreSamplesColsAreLocations = pe.getSnpData(...) %>%
      unite(alleles, c(allele1, allele2), sep = "") %>%
      select(-one_of("avgGenCallScore")) %>%
      spread(illuminaSnpLocation, alleles)

    withoutNonCalls = pe.noncalls.dropEntireColumn(rowsAreSamplesColsAreLocations)
    
    theMatrix = apply(withoutNonCalls[,-1], 2, function (col) {
      tb = table(col)
      dominant = names(which.max(tb))
      recessive = names(which.min(tb))
      if (dominant == recessive | !(substr(recessive, 1, 1) == substr(recessive, 2, 2)))
        recessive = "Z" # This just ensures we only set recessive if you have both
      as.integer(gsub("[^\\d\\-]{2}", 0, gsub(dominant, 1, gsub(recessive, -1, col, fixed = TRUE), fixed = TRUE)))
    })
    
    theMatrix %>% pe.dropMonomorphicColumns
  })
}

pe.getIdsForPatientsWhoQualified = function () cache.get("clinicalIds", function () pe.getSnpData()$sampleId %>% unique)

pe.getClinicalInfo = function () {
  clinicalIds <- pe.getIdsForPatientsWhoQualified()
  read_csv(pe.clinicalDataPath, col_types = cols()) %>%
    filter(patientId %in% clinicalIds) %>%
    arrange(patientId) %>%
    mutate(mixed = as.integer(!is.na(ethnicity2)), ethnicityCategory = as.factor(ethnicityCategory))
}

pe.getSnpsWithChipInfo = function () pe.getSnpData() %>% cbind(pe.getMetaboChip() %>% select(-one_of(c("SNP", "Name", "Ploidy", "Species"))))

pe.getMetaboChip = function () {
  cache.get("metaboChip", function () {
    read_csv(pe.illuminaPath, skip=7, n_max = 196725, col_types = cols(Chr = col_character())) %>% mutate(Chr = as.factor(Chr))
  })
}

pe.getClinicalAsNumerics = function () {
  pe.getClinicalInfo() %>%
    select(preeclampsia, momAge, gestAge, ethnicityCategory, primaGravida, gdm, iugr, abruption, daughter, babyWgtGram ) %>%
    mutate(ethnicityCategory = as.integer(ethnicityCategory))
}

pe.getDuplicateMeasurementLocations = function () {
  metaboChip = pe.getMetaboChip()
  allPositions = paste(metaboChip$Chr, metaboChip$MapInfo)
  
  metaboChip[duplicated(allPositions) | duplicated(allPositions, fromLast = T),] %>%
    select(Chr, MapInfo) %>% mutate(id = paste(Chr, MapInfo)) %>% arrange(Chr, MapInfo)
}

pe.getSimpleSNPFormatForPhyloTree = function() {
  cache.get("simpleSNP", function () {
    "%notin%" <- Negate("%in%")

    fullSet = pe.getSnpsWithChipInfo()
    duplicates = pe.getDuplicateMeasurementLocations()
    simpleSNP = fullSet %>%
      select(Chr, MapInfo, sampleId, allele1) %>%
      mutate(id = paste(Chr, MapInfo)) %>% filter(id %notin% duplicates$id) %>% select(-one_of("id")) %>%
      filter(!grepl('-', allele1))
    
    colnames(simpleSNP) = c("#Chrom", "Pos", "Sample", "Allele")
    simpleSNP = simpleSNP %>% spread(Sample, Allele, sep = "") %>% na.omit
    simpleSNP$`#Chrom` = as.integer(pe.recodeChrom(simpleSNP$`#Chrom`))
    simpleSNP
  })
}

pe.exploreLoci = function (loci, dropNonCalls = FALSE) {
  res = pe.getSnpData() %>% filter(illuminaSnpLocation == loci) %>% cbind(pe.getClinicalInfo()) %>%
    unite(alleles, c(allele1, allele2), sep = "") %>% group_by(alleles, preeclampsia) %>%
    summarise(n = n()) %>% arrange(preeclampsia) %>% spread(alleles, n) %>% as.data.frame
  rownames(res) = c("control", "eoPreclampsia")
  res$preeclampsia = NULL
  if (dropNonCalls & ncol(res) == 4)
    res = res[,-1]
  res[is.na(res)] <- 0
  res
}

pe.getSnpTable = function (patients = pe.getClinicalInfo(), ...) {
  cache.get("snpTable", function () {
    samps = pe.getSamplesAsRows(...) %>% filter(sampleId %in% patients$patientId)
    
    samplesAsRows = patients %>% select(preeclampsia) %>% cbind(samps) %>%
      pe.dropColumnsThatContainNA %>% pe.dropMonomorphicColumns
    
    setupSNP(samplesAsRows, 3:ncol(samplesAsRows),sep="")
  })
}

pe.getSnpTableFilip = function (...) pe.getSnpTable(pe.getClinicalInfo() %>% filter(ethnicityCategory == "Filipino"), ...)

pe.makeFilesForSnpQCprogram = function (outputDir = "snpqc/data") {
  # File 1
  clinical = pe.getClinicalInfo()
  clinical %>% mutate(ID = patientId) %>% select(ID, preeclampsia) %>% write_tsv(glue("{outputDir}/sampleMap.txt"))
  
  # File 2
  snpdata =  pe.getSnpData()
  metabo = pe.getMetaboChip() %>% arrange(Name)
  metabo$Name = gsub("[\\:\\-\\_]", "\\:", metabo$Name)
  snpdata$illuminaSnpLocation = gsub("[\\:\\-\\_]", "\\:", snpdata$illuminaSnpLocation)
  metabo %>% select(Chr, MapInfo, Name) %>%
    mutate(Chromosome = Chr, Position = MapInfo) %>%
    select(Name,Chromosome,Position) %>% write_tsv(glue("{outputDir}/SNPmap.txt"))

  # File 3: Needs SNPs in A/B format like so:
  # SNP Name	Sample ID	Allele1 - AB	Allele2 - AB GC Score
  # snp1	sample1	A	B	0.8446
  # snp2	sample1	B	B	0.9629
  # snp3	sample1	B	B	0.9484
  topAlleles = snpdata %>% mutate(allele = allele1) %>% select(illuminaSnpLocation, allele) %>%
    rbind(snpdata %>% mutate(allele = allele1) %>% select(illuminaSnpLocation, allele)) %>%
    group_by(illuminaSnpLocation) %>% count(allele) %>% top_n(2)
  
  topAllelesCalled = topAlleles %>% filter(allele != "-") %>% arrange(desc(n))
  topAllelesCalled = topAllelesCalled[!duplicated(topAllelesCalled$illuminaSnpLocation),]
  
  # switch to A/B
  snps = left_join(snpdata, topAllelesCalled, by = "illuminaSnpLocation")
  snps$alleleA = ifelse(snps$allele1 == "-", "-", ifelse(snps$allele1 == snps$allele, "A", "B"))
  snps$alleleB = ifelse(snps$allele2 == "-", "-", ifelse(snps$allele2 == snps$allele, "A", "B"))
  snps %>% select(illuminaSnpLocation, sampleId, alleleA, alleleB, avgGenCallScore) %>% write_tsv(glue("{outputDir}/SNPsample.txt"))
}

pe.recodeChrom = function (chr) recode(chr, "X" = "23", "XY" = "24", "Y" = "25", "Mt" = "26")

pe.getPValsAll = function (suffix = "") {
  cache.get(glue("pvalues{suffix}"), function () {
    getPvals = function (key) {
      cache.get(key, function () stop(glue("{key} error: You need to use thread.r to compute pvals with one model per thread, otherwise will take days.")))
    }
    
    assocResults = cbind(getPvals(glue("dominant{suffix}")),
                         getPvals(glue("codominant{suffix}")),
                         getPvals(glue("recessive{suffix}")),
                         getPvals(glue("log-additive{suffix}")),
                         getPvals(glue("overdominant{suffix}")))[,c(2,4,6,8,10)]

    pvalueResults = assocResults
    pvalueResults$minPval <- as.numeric(apply(pvalueResults, 1, function (row) {
      values = na.omit(row[2:6])
      if (length(values) == 0)
        return(1)
      min(values)
    }))
    
    pvalueResults$bonferroni = p.adjust(pvalueResults$minPval, method = "bonferroni")
    pvalueResults$holm = p.adjust(pvalueResults$minPval, method = "holm")

    pvalueResults$location = gsub("[\\:\\-\\_]", "\\.", rownames(pvalueResults))
    metabo = pe.getMetaboChip() %>% arrange(Name)
    metabo$location = gsub("[\\:\\-\\_]", "\\.", metabo$Name)
    joined = left_join(pvalueResults, metabo %>% select(location, Chr), by = "location")

    stopifnot(sum(is.na(joined$Chr)) == 0)
    
    pvalueResults$chr = as.integer(pe.recodeChrom(as.character(joined$Chr)))
    pvalueResults
  })
}
