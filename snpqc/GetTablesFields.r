require("RSQLite")

TableNames=dbListTables(dbcon) # get tables - not much point now but might want to make more generic in the future
tSamples="samples" #TableNames[Sindex] # samples
tSNPmap="snpmap" #TableNames[SMindex]  # snpmap
tSNP="snps"  #TableNames[SNPindex]     # snps

#get names used in columns in different tables
# sample info
FieldNames=dbListFields(dbcon,tSamples) # fields in table
fIDS=FieldNames[sampleindexR]

# map info
FieldNames=dbListFields(dbcon,tSNPmap) # fields in table
fsnpM=FieldNames[snpindexR]
fchromM=FieldNames[chromindex]
fposM=FieldNames[bpindex]

# genotypes info
FieldNames=dbListFields(dbcon,tSNP) # fields in table
fID=FieldNames[genosampleindex]
fsnp=FieldNames[genosnpindex]
fall1=FieldNames[genoallele1index]
fall2=FieldNames[genoallele2index]

fgc=fall2 # if no gcscores available, simply uses the allele2 index as a placeholder and then uses a value of One for the gcscore (largely beats the point, but does some population parameter filtering).
hasGC=F
if (genogcindex > 0)
{
  fgc=FieldNames[genogcindex]
  hasGC=T
}

# not implemented
fx=0
fy=0
if (genoxindex >0) fx=FieldNames[genoxindex]
if (genoyindex >0) fy=FieldNames[genoyindex]

hasXY=F
if (genoxindex>0 & genoyindex>0) hasXY=T

rm(TableNames,FieldNames,sampleindexR,snpindexR,chromindex,bpindex,genosampleindex,genosnpindex,genoallele1index,genoallele2index,genoxindex,genoyindex)

