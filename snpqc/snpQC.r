#### Entry point for snp QC script

require("gplots")
require("gtools") 
require("lattice")
require("RSQLite")
require("compiler")
require("snowfall")
require("gdata")

enableJIT(3) 

.maxcores=32 # increase max number of cores
.cpu=8 # setup number of cores in cluster
.machineName="localhost"

#separator="\t" # flat files tab separated - could make more flexible in future if needed.
# numHits=10000 # number of snps per iteration 

cat("\n\nsnpQC ready to run\ncall runQC(runParam='path/name of parameters file',numHits=10000,DBin=TRUE,qc=TRUE,filter=TRUE,report=TRUE,DBout=TRUE)\n") # start time

source("QCfunctions.r")  # load QC functions
source("SummaryAndPlots.r") # load plotting and summary function

useYates=T
buildGRM=F
Text=""
separator=""
numHits=10000
.Data=NULL

runQC=function(runParam="RunParameters.txt",numHits=10000,DBin=TRUE,qc=TRUE,filter=TRUE,report=TRUE,DBout=TRUE,useYates=TRUE)
{               
  options(warn=-1)
  cat(paste("job started on",date(),"\n")) # start time
  Text<<-.readParams(runParam)
  numHits<<-numHits
  separator<<-"\t"
  source("ReadParameters.r")
  if (DBin==TRUE) .runDB(runParam)
  if (qc==TRUE) .runQC()
  if (filter==TRUE) .runPlotsAndFilters()
  if (report==TRUE) .runReport()
  if (DBout==TRUE).runDBresults()
  if (buildGRM==T) source("BuildGRM.r") # can be called as a pass through with all set to FALSE 
  cat(paste("job completed on",date(),"\n")) # end time
  source("Cleanup.r") 
  options(warn=0)  
}

.readParams=function(runParam)
{
  Con=file(paste(runParam,sep=""),"r")
  Text=readLines(Con)
  close(Con)   
  return (Text)
}

.runDB=function(runParam) # builds DB function
{
  source("ReadParameters.r")
  source("BuildDB.r")    
}

.runQC=function() # runs map, snp and sample qc
{
  source("ReadParameters.r")
  source("RunQC.r")   
}

.runPlotsAndFilters=function() # build summaries and plots - only sample summary not generated here
{
  source("ReadParameters.r")   
  .summaryAndPlots()    
}

.runReport=function() # build Tex report
{
  source("ReadParameters.r")
  Sweave("SNPreport.rnw")
  
  # see if can run/find latex
  system2(command="pdflatex", args="SNPreport.tex", wait=TRUE,stdout=FALSE)  
  system2(command="pdflatex", args="SNPreport.tex", wait=TRUE, stdout=FALSE)
  system2(command="pdflatex", args="SNPreport.tex", wait=TRUE,stdout=FALSE)
    
  if (file.exists("SNPreport.log")==TRUE) file.remove("SNPreport.log")
  if (file.exists("SNPreport.aux")==TRUE)file.remove("SNPreport.aux")
  
  if (file.exists("SNPreport.tex")==TRUE) 
  {
    file.copy("SNPreport.tex",to=paste(reportname,"snpQCreport.tex",sep="/"),overwrite = TRUE)
    file.remove("SNPreport.tex")    
  }
  
  if (file.exists("SNPreport.pdf")==TRUE) 
  {
    file.copy("SNPreport.pdf",to=paste(reportname,"snpQCreport.pdf",sep="/"),overwrite = TRUE)
    file.remove("SNPreport.pdf")    
  }  
}

.runDBresults=function() # populate results DB
{
  source("ReadParameters.r")
  source("BuildResultsDB.r")
}  




