require("RSQLite")
dbcon=dbConnect(dbDriver("SQLite"), dbname=dbname) # connect to database

# crank up speed - might still not be optimal
dbGetQuery(dbcon, "pragma cache_size = 800000")
dbGetQuery(dbcon, "pragma synchronous = 0")
dbGetQuery(dbcon, "pragma temp_store = 2")
dbGetQuery(dbcon, "pragma journal_mode = off")
dbGetQuery(dbcon, "pragma page_size=65536")
dbGetQuery(dbcon, "pragma locking_mode=exclusive")

source("GetTablesFields.r") # get indices of data in DB - based on user info

# get samples from DB
animids=as.vector(dbGetQuery(dbcon,paste("select distinct ",fIDS," from ",tSamples,sep=""))[,1]) # get all sample ids
animids=mixedsort(animids)
numanim=length(animids) # number of animals

snps=.mapStats() # mapping stats
numsnp=length(snps) # number of snps

.snpStats() # snp stats
numgeno=.sampleStats() # sample stats
dbDisconnect(dbcon) # close DB connection
.sampleCorrelation() # sample correlations