require("gplots")

# index for parameter file info - make these more robust in future (XML or string match)

# flat files with raw data
snpMapPar=2
samplePar=3
genoPar=4

# column ids of data (flat files and DB)
genosnpindexPar=6				 # index of SNP in genotype file - must match names in SNP map file or database column
genosampleindexPar=7		 # index of sample identifiers in genotype file - must match names in sample map file or database column
genoallele1indexPar=8	 # index of SNP allele1 in genotype file - currently only uses A/B allele calls or database column
genoallele2indexPar=9	 # index of SNP allele2 in genotype file - currently only uses A/B allele calls or database column
genogcindexPar=10				 # index of gcscore in genotype file or database column
sampleindexPar=11				 # index of sample identifiers in sample map file or database column
snpindexPar=12				   # index of SNP names in SNP map file or database column
chromindexPar=13				 # index of chromosome in SNP map file or database column
bpindexPar=14				     # index of base pair position in SNP map file or database column

# will need these for the pictures
genoxindexPar=45
genoyindexPar=46

# number of lines of header to skip when building DB
skipgenoPar=15
skipsamplePar=16
skipsnpPar=17

# missing symbol for alleles
misPar=18

# filtering criteria for SNP
snpPar=20:28

# filtering criteria for SNP
sampleFilterPar=29:31

# snps to exclude from mapping info - user defined
mapsnpPar=33

# parameters for output files
outPar=35:39

# database info
dbPar=41 # db name
dbDirFlagPar=42 # use default location yes/no
dbDir=43 # custom location for db

# info for report - project name and authors
projPar=48
authorPar=49

##################################################################
# read in parameters file
# get database path and name - create directories for output
dbname=trim(gsub("\t","",unlist(strsplit(Text[dbPar],"#"))[1])) # get rid of common junk - tabs and empty spaces
dbLocBool=toupper(trim(gsub("\t","",unlist(strsplit(Text[dbDirFlagPar],"#"))[1])))
dbDirPath=""
if (dbLocBool == FALSE)
{
  dbDirPath=trim(gsub("\t","",unlist(strsplit(Text[dbDir],"#"))[1]))
  if (file.exists(dbDirPath)==FALSE)
  {
    cat("The path for the database does not exist. Creating...\n\n")
    dir.create(dbDirPath)
  }
} else
{
  dbDirPath="databases"
  if (file.exists("databases")==FALSE) dir.create("databases")
}
outputname=paste(dbDirPath,"/",dbname,"_output",sep="")
if (file.exists(outputname)==FALSE) dir.create(outputname)
reportname=paste(dbDirPath,"/",dbname,"_report",sep="")
if (file.exists(reportname)==FALSE) dir.create(reportname)
dbname=paste(dbDirPath,"/",dbname,".db",sep="")

outputname=gsub("\\\\","/",outputname)  # replace backslash for LateX
reportname=gsub("\\\\","/",reportname)

# get raw data file path and name
snpPath=trim(gsub("\t","",unlist(strsplit(Text[snpMapPar],"#"))[1]))
samplePath=trim(gsub("\t","",unlist(strsplit(Text[samplePar],"#"))[1]))
genoPath=trim(gsub("\t","",unlist(strsplit(Text[genoPar],"#"))[1]))

#   number of header lines to skip in raw data
skipgeno=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[skipgenoPar],"#"))[1])))
skipsample=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[skipsamplePar],"#"))[1])))
skipsnp=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[skipsnpPar],"#"))[1])))

genosnpindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genosnpindexPar],"#"))[1])))
genosampleindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genosampleindexPar],"#"))[1])))
genoallele1index=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genoallele1indexPar],"#"))[1])))
genoallele2index=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genoallele2indexPar],"#"))[1])))
genogcindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genogcindexPar],"#"))[1])))
sampleindexR=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[sampleindexPar],"#"))[1])))
snpindexR=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpindexPar],"#"))[1])))
chromindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[chromindexPar],"#"))[1])))
bpindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[bpindexPar],"#"))[1])))

# X/Y just for figure - good/bad
genoxindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genoxindexPar],"#"))[1])))
genoyindex=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[genoyindexPar],"#"))[1])))

 # parameters for SNP filtering
# could be on a single line - but easier to see comments
snparams=numeric(9)
snparams[1]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[1]],"#"))[1])))  # over X percent genotyping fail - snp criterion
snparams[2]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[2]],"#"))[1])))  # median call rates smaller than X  - snp criterion
snparams[3]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[3]],"#"))[1])))  # all GC scores zero or NA (do not change) - snp criterion
snparams[4]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[4]],"#"))[1])))  # GC < X in Y (below) percent of samples - snp criterion
snparams[5]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[5]],"#"))[1])))  # GC < X (above) in Y percent of samples - snp criterion
snparams[6]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[6]],"#"))[1])))  # MAF zero - snps fixed across all samples, use negative value not to exclude snps that are not segregating
snparams[7]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[7]],"#"))[1])))  # minor allele frequency MAF < X - snp criterion
snparams[8]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[8]],"#"))[1])))  # heterozygosity deviation in number of standard deviations X - snp criterion
snparams[9]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[snpPar[9]],"#"))[1])))  # Bonferroni corrected Hardy-Weinberg threshold X - snp criterion

GCcutoff=snparams[4] # cutoff snps with gc scores under this value

# missing data symbol - R version uses only the first one
missingSymbol=trim(unlist(strsplit(Text[misPar],"#")))[1]
missingSymbol=unique(trim(unlist(strsplit(missingSymbol,"\t"))))
missingSymbol=missingSymbol[1]

# snps to exclude from map info
mapexcludeSNP=trim(unlist(strsplit(Text[mapsnpPar],"#")))[1]
mapexcludeSNP=unique(trim(unlist(strsplit(mapexcludeSNP,"\t"))))
mapexcludeSNP=mapexcludeSNP[which(mapexcludeSNP!="")]

# output parameters
excludeReshaped=toupper(trim(gsub("\t","",unlist(strsplit(Text[outPar[1]],"#"))[1]))) #exclude all rejected SNP and samples from 'reshaped.txt' file - T/F
setReshaped9=toupper(trim(gsub("\t","",unlist(strsplit(Text[outPar[2]],"#"))[1]))) # set all rejected samples and SNP to missing value
missingReshaped=trim(gsub("\t","",unlist(strsplit(Text[outPar[3]],"#"))[1])) #symbol for missing/rejected values

genoReshaped=trim(unlist(strsplit(Text[outPar[4]],"#")))[1]
genoReshaped=unique(trim(unlist(strsplit(genoReshaped,"\t"))))
genoReshaped=genoReshaped[which(genoReshaped!="")] #values for genotypes in reshaped file

buildGRM=toupper(trim(gsub("\t","",unlist(strsplit(Text[outPar[5]],"#"))[1]))) # build GRM matrix T/F

# parameters for removal of bad samples
sampleparams=numeric(3)
sampleparams[1]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[sampleFilterPar[1]],"#"))[1])))  # call rates lower than X - sample criterion
sampleparams[2]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[sampleFilterPar[2]],"#"))[1])))  # eterozygosity deviation in number of standard deviations X - sample criterion
sampleparams[3]=as.numeric(trim(gsub("\t","",unlist(strsplit(Text[sampleFilterPar[3]],"#"))[1])))  # correlation between samples (just shown but not used to exclude samples)

# info for report - project name and authors
projectname=trim(unlist(strsplit(Text[projPar],"#")))[1]
authors=trim(unlist(strsplit(Text[authorPar],"#")))[1]


# clean up
rm(dbLocBool,dbDirPath, snpMapPar,samplePar,genoPar, genosnpindexPar, genosampleindexPar,genoallele1indexPar,
    genoallele2indexPar,genogcindexPar,sampleindexPar, snpindexPar, chromindexPar,bpindexPar, genoxindexPar,
    genoyindexPar,skipgenoPar,skipsamplePar,skipsnpPar, misPar,snpPar, sampleFilterPar, mapsnpPar,
    outPar,dbPar,dbDirFlagPar,dbDir, projPar, authorPar)

