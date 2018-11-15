####### start building DB ########################
require("RSQLite")

begin=Sys.time() # keep track of how long it takes
cat("building DB...\n\n")

# create and connect to new database
dbcon=dbConnect(dbDriver("SQLite"), dbname = dbname)

# speed up DB - more risky
dbGetQuery(dbcon, "pragma cache_size = 200000")
dbGetQuery(dbcon, "pragma synchronous = 0")
dbGetQuery(dbcon, "pragma temp_store = 2")
dbGetQuery(dbcon, "pragma journal_mode = off")
dbGetQuery(dbcon, "pragma page_size=65536")
dbGetQuery(dbcon, "pragma locking_mode=exclusive")

#EOL="\n"
#if (Sys.info()[1]=="Windows") EOL="\r\n" # avoid problems with return across OS

######## # populate tables  ##################
# test each file for EOL - has caused problems with mishing and mashing of files
EOL="\n"
ret=readChar(samplePath,5000)
if (length(grep("\r\n",ret))>0) EOL="\r\n"
dbWriteTable(dbcon,"samples",samplePath,append=T,header=F,skip=skipsample,sep=separator,eol=EOL) # change skip to the number of headers before the actual data

EOL="\n"
ret=readChar(snpPath,5000)
if (length(grep("\r\n",ret))>0) EOL="\r\n"
dbWriteTable(dbcon,"snpmap",snpPath,overwrite=T,header=F,skip=skipsnp,sep=separator,eol=EOL) # change skip to the number of headers before the actual data

EOL="\n"
ret=readChar(genoPath,5000)
if (length(grep("\r\n",ret))>0) EOL="\r\n"
dbWriteTable(dbcon,"snps",genoPath,append=T,header=F,skip=skipgeno,sep=separator,eol=EOL) # change skip to the number of headers before the actual data
##########################################################################################
rm(snpPath,samplePath,genoPath,skipgeno,skipsample,skipsnp) # clean up

print(Sys.time()-begin) # show runtime

####### start building DB indices########################
begin=Sys.time() # keep track of how long it takes
dbGetQuery(dbcon, "pragma cache_size = 2000") # apparently smaller is better for indexing - need to test

source("GetTablesFields.r")

cat("building indices...\n\n")

cat("indexing samples...\n\n")

hold=dbGetQuery(dbcon,"select * from sqlite_master where type='index'") # get all indexes

if(length(which(hold$name=="sample_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX sample_idx")) # remove old index
dbGetQuery(dbcon, paste("CREATE INDEX sample_idx ON samples(",fIDS,")",sep="")) # index samples

cat("indexing snp map...\n\n")
if(length(which(hold$name=="chromosome_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX chromosome_idx")) # remove old index
if(length(which(hold$name=="snpmap_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX snpmap_idx")) # remove old index
dbGetQuery(dbcon, paste("CREATE INDEX chromosome_idx ON snpmap(",fchromM,")",sep="")) # index chromosome
dbGetQuery(dbcon,paste("CREATE INDEX snpmap_idx ON snpmap(",fsnpM,")",sep="")) # index snp

cat("indexing genotypes (can take a very long time)...\n\n")
if(length(which(hold$name=="snp_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX snp_idx")) # remove old index
if(length(which(hold$name=="ID_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX ID_idx")) # remove old index
dbGetQuery(dbcon,paste("CREATE INDEX snp_idx ON snps(",fID,")",sep=""))  # index snp on genotypes
dbGetQuery(dbcon,paste("CREATE INDEX ID_idx ON snps(",fsnp,")",sep=""))  # index samples on genotypes

cat("indexing complete!!\n\n")

dbDisconnect(dbcon) # close DB connection
print(Sys.time()-begin) # show runtime


