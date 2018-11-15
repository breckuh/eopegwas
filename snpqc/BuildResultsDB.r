####### start building DB ########################
require("RSQLite")

begin=Sys.time() # keep track of how long it takes
cat("adding results to DB...\n\n")

# create and connect to new database
dbcon=dbConnect(dbDriver("SQLite"), dbname = dbname)

# speed up DB - more risky
dbGetQuery(dbcon, "pragma cache_size = 200000")
dbGetQuery(dbcon, "pragma synchronous = 0")
dbGetQuery(dbcon, "pragma temp_store = 2")
dbGetQuery(dbcon, "pragma journal_mode = off")
dbGetQuery(dbcon, "pragma page_size=65536")
dbGetQuery(dbcon, "pragma locking_mode=exclusive")

EOL="\n"
if (Sys.info()[1]=="Windows") EOL="\r\n" # avoid problems with return across OS

######## # populate tables  ##################
dbWriteTable(dbcon,"snpfilter",paste(outputname,"snpfilter.txt",sep="/"),overwrite=T,skip=0,header=TRUE,sep=separator,eol=EOL)
dbWriteTable(dbcon,"samplefilter",paste(outputname,"samplefilter.txt",sep="/"),overwrite=T,skip=0,header=TRUE,sep=separator,eol=EOL)
##########################################################################################

print(Sys.time()-begin) # show runtime

####### start building DB indices########################
begin=Sys.time() # keep track of how long it takes
dbGetQuery(dbcon, "pragma cache_size = 2000") # apparently smaller is better for indexing - need to test

cat("building indices...\n\n")

cat("indexing samples...\n\n")

hold=dbGetQuery(dbcon,"select * from sqlite_master where type='index'") # get all indexes

if(length(which(hold$name=="sampleQC_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX sampleQC_idx")) # remove old index  
dbGetQuery(dbcon, "CREATE INDEX sampleQC_idx ON samplefilter(good)") # index good/bad sample

cat("indexing snps...\n\n")
if(length(which(hold$name=="snpQC_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX snpQC_idx")) # remove old index  
if(length(which(hold$name=="snpMapped_idx"))>0) dbGetQuery(dbcon, paste("DROP INDEX snpMapped_idx")) # remove old index  
dbGetQuery(dbcon, "CREATE INDEX snpQC_idx ON snpfilter(good)") # index good/bad SNP
dbGetQuery(dbcon, "CREATE INDEX snpMapped_idx ON snpfilter(mapped)") # index mapping exclusion for SNP

cat("indexing complete!!\n\n")

dbDisconnect(dbcon) # close DB connection
rm(hold)

print(Sys.time()-begin) # show runtime


