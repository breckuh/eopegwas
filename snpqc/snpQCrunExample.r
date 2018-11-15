# example using R version

# change working directory to where the snpQC scripts are on disk
# setwd("") # e.g. path with snpQC scripts
source("snpQC.r") # load scripts
runParam="RunParam.txt" # path and name of the parameters file (see example in dataset)
runQC(runParam) # runs the full QC
#runQC(runParam,numHits=50000,DBin=TRUE,qc=TRUE,filter=TRUE,report=TRUE,DBout=TRUE) # all parameters in the function
# that's all


###### some details ########
# user's can define what parts of the QC to run
# the full list of parameters is
#runQC(runParam,numHits=50000,DBin=TRUE,qc=TRUE,filter=TRUE,report=TRUE,DBout=TRUE)

# runParam = path and name of the parameters file
# numHits = how many SNP to read in each time, more is faster but more demanding on the computational resources
# DBin = build database of genotypes
# qc = run QC metrics
# filter = summarize results and build plots for report
# report = make PDF report
# DBout = write QC results to database

# note: none of DBin, qc, report or DBout can be run if the previous steps have not been run before (as in never before).
# you'll need the database to run the qc
# you'll need qc results to build the summaries/plots
# you'll need summaries/plots to write the report
# you can actually skip filter and report and go straight from the qc to writing the results to the database

# hidden parameters that might be useful
# .maxcores=32 # increase max number of cores
# .cpu=8 # setup number of cores in cluster, can use more or less - if more than 32 need to increase .maxcores
# .machineName="localhost" # change to run on a remote machine, will probably need to sort out ssh tunneling to get it to work
# useYates=TRUE # True or false for using Yates correction on the chi-square values

# the program is pretty slow - particularly to build the DB, but should run on just about anything 
# a fast version in C# is also available and can bypass the DB is this is not needed. Usage is exactly the same. Will only run in Windows 64 bits with .Net 4.5



###################################################################################################
#################### installing required packages ################################################
#################################################################################################

# snpQC uses a few packages - you can run the following line to install them:
#install.packages(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"),dependencies=T)
# and then choose a mirror to download from (a list should pop up)

# if you are running R without a graphical interface, use: 
#setRepositories(graphics=F, ind=1:2) # set repositories for CRAN
#local({r = getOption("repos"); r["CRAN"] = "http://cran.ms.unimelb.edu.au"; options(repos=r)}) # set to Melbourne mirror (or any mirror you prefer)
#install.packages(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"),dependencies=T)

# if you are running on a server you might have problems with permissions to install packages - create a personal library, install the packages there and set the path to it 
# to install:
#myRlib="/home/user/Rlib" # for example
#install.packages(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"),dependencies=T,lib="myRlib")

# and to use the packages, set the path with 
#.libPaths("/home/user/Rlib/") # before running snpQC  
