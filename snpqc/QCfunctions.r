# SNP analysis - preprocessing
# for illumina data
require("compiler")
enableJIT(3) 

############################
# Mapping statistics ######
##########################
.mapStats=function()
{
  # get snps and sort by chromosome X bp
  begin=Sys.time() # keep track of how long it takes 
  cat("building ordered map...\n")
  
  require("gtools") # for just natural number ordering of chromosomes
      
  map=dbGetQuery(dbcon,paste("select ",fsnpM,", cast (",fchromM," as char),",fposM," from ",tSNPmap,sep=""))
  
  map=data.frame(snp=map[,1],chrom=map[,2],loc=map[,3])
  map[,1]=as.character(map[,1])
  map[,2]=as.factor(map[,2])
  map[,3]=as.numeric(as.character(map[,3]))
  map=map[mixedorder(map[,3]),] # sort BP
  map=map[mixedorder(map[,2]),] # sort chromosomes
 
  snps=as.character(map[,1])
  write.table(map,paste(outputname,"orderedmap.txt",sep="/"),quote=F,sep="\t",row.names=F)
  rm(map) # clean up
  
  cat("map ordering complete.\n")
  print(Sys.time()-begin) # show runtime
  return(snps)
} # end map function



#########################
#### per SNP QC ########
#######################
.snpStats=function()
{
  # numHits=10000 # number of snps per round
  begin=Sys.time() # keep track of how long it takes 
  cat("running SNP filtering...\n")
  
  require("snowfall")
  #machineName="localhost" # values from RunSNPAnalysis.r - can add to parameter file in future
  #cpu=8
  #maxcores=32
  sfSetMaxCPUs(number=.maxcores) # increase max number of cores
  sfInit(parallel=TRUE,cpus=.cpu,socketHosts=rep(.machineName,.cpu)) # setup number of cores in cluster
    
  SummarizeGCinfo=function(tempgc)
  {
    gchold=numeric(12)
    gchold[1]=length(tempgc[which(tempgc==0)]) # number of gc score=0
    gchold[2]=length(tempgc)-as.numeric(as.character(gchold[1])) # number of not zeros    
    gchold[3]=min(tempgc) # minimum
    gchold[4]=max(tempgc) # maximum
    gchold[5]=mean(tempgc) # mean
    gchold[6]=median(tempgc) # median
    gchold[7]=length(tempgc[tempgc<0.1]) # <0.1
    gchold[8]=length(tempgc[tempgc>=0.1 & tempgc<0.5]) # <0.5
    gchold[9]=length(tempgc[tempgc>=0.5 & tempgc<0.6]) # <0.6
    gchold[10]=length(tempgc[tempgc>=0.6 & tempgc<0.9]) # <0.9
    gchold[11]=length(tempgc[tempgc>=0.9]) # >=0.9
    gchold[12]=length(tempgc[tempgc<snparams[4]]) # < user defined value 
    return (gchold)
  } # end GC function
  
  # get data in blocks of snps   
  firstSNP=1 # counter 1
  lastSNP=numHits # counter 2
  
  # initialize text files 
  Text=paste("Missing","Called","Min","Max","Mean","Median","0.0_0.1","0.1_0.5","0.5_0.6","0.6_0.9","0.9_1",paste("0_",snparams[4],sep=""),sep="\t")
  write(Text,paste(outputname,"gcstats.txt",sep="/"),append=F)
  
  Text=paste("MissAllele","A","B","MissGeno","AA","AB","BB",sep="\t")
  write(Text,paste(outputname,"snpfreqs.txt",sep="/"),append=F)
  write(Text,paste(outputname,"snpcounts.txt",sep="/"),append=F)
  
  Text=paste("chisquare","pvalue",sep="\t")
  write(Text,paste(outputname,"hwstats.txt",sep="/"),append=F)
  
  rm(Text) # clean up
  
  # start SNP loop
  numIter=floor(numsnp/numHits)
  lastOnes=numsnp-numIter*numHits
  if (lastOnes>0) numIter=numIter+1
  
  for (i in 1:numIter)
  { 
    if(i==numIter) lastSNP=numsnp # last loop get just the missing ones
    cat(paste("evaluating SNPS:",firstSNP,"-",lastSNP,"\n"))
    # get data
    getSNP=snps[firstSNP:lastSNP]
    curNumSNP=length(getSNP)
    dbSNPQuery=paste(getSNP,collapse="','")
    hold=dbGetQuery(dbcon,paste("select ",paste(fID,fsnp,fall1,fall2,fgc,sep=",")," from ",tSNP," where ",fsnp," in ('",dbSNPQuery,"')",sep=""))
    names(hold)=c("ID","snp","allele1","allele2","gcscore")
    
    idindex=hold$snp[seq(from=1,to=curNumSNP*numanim,by=numanim)]  # fix up order - DB does not return data in order of snps requested
    idindex=match(getSNP,idindex) # match with requested snps
      
    if (hasGC==FALSE) hold$gcscore=1 # set gcscore to one for all   
    
    # convert to numeric from AB format - will not work with other formats
    # add something here to fix up in future - maybe based on the map info, but needs to recode on a per snp basis 
    geno=paste(hold$allele1,hold$allele2,sep="")
    geno[geno=="AA"]=0
    geno[geno=="AB"]=1
    geno[geno=="BA"]=1
    geno[geno=="BB"]=2
    geno[geno==paste(missingSymbol,missingSymbol,sep="")]=NA
    geno=matrix(as.numeric(geno),curNumSNP,numanim,byrow=T) # numeric matrix of genotypes
        
    # gcscore summary
    gcscore=matrix(as.numeric(as.character(hold$gcscore)),curNumSNP,numanim,byrow=T) # make sure numeric - missing might be something different from NA, will convert to NA but give warning
    gcscore[is.na(gcscore)==T]=0 # set NAs to zero
    gcscore[is.nan(gcscore)==T]=0 # set NANs to zero
    sfExport(list=list("snparams"))
    gcstats=matrix(t(sfApply(gcscore,1,SummarizeGCinfo)),curNumSNP,12) 
    gcstats=gcstats[idindex,]
    write(t(cbind(snps[firstSNP:lastSNP],gcstats)),paste(outputname,"gcstats.txt",sep="/"),ncolumns=13,sep="\t",append=T) # write output
    
    if (hasGC==TRUE) geno[gcscore<GCcutoff]=NA  # set genotypes below threshold as missing
    # genotypes and alleles - counts and frequencies
    # genotypes
    misgeno=unlist(sfApply(geno,1,function (x) length(which(is.na(x)==T)))) # number of missing genotypes
    genotot=(numanim-misgeno) # total number of called genotypes 
    AAcount=unlist(sfApply(geno,1,function (x) length(which(x==0)))) # number of AA genotypes
    ABcount=unlist(sfApply(geno,1,function (x) length(which(x==1)))) # number of AB genotypes
    BBcount=genotot-(AAcount+ABcount) # number of BB genotypes
    #alleles
    alleletot=genotot*2 # total number of called alleles
    Bcount=rowSums(geno,na.rm=T) # B count
    Acount=alleletot-Bcount # A count
    q=rowSums(geno,na.rm=T)/alleletot # B frequency
    p=1-q  # A frequency 
    
    snpcounts=matrix(c(misgeno*2,Acount,Bcount,misgeno,AAcount,ABcount,BBcount),curNumSNP,7)
    snpcounts=snpcounts[idindex,]
    write(t(cbind(snps[firstSNP:lastSNP],snpcounts)),paste(outputname,"snpcounts.txt",sep="/"),ncolumns=8,sep="\t",append=T) # write output
    
    snpfreqs=matrix(c((misgeno*2)/(numanim*2),p,q,misgeno/numanim,AAcount/genotot,ABcount/genotot,BBcount/genotot),curNumSNP,7)
    snpfreqs=snpfreqs[idindex,]
    write(t(cbind(snps[firstSNP:lastSNP],snpfreqs)),paste(outputname,"snpfreqs.txt",sep="/"),ncolumns=8,sep="\t",append=T)  # write output
    
    # HW equilibrium 
    Exp=matrix(c(p^2,2*p*q,q^2),curNumSNP,3)
    Exp=Exp*genotot
    Exp=t(Exp) 
    
    Obs=matrix(c(AAcount,ABcount,BBcount),3,curNumSNP,byrow=T)
    xtot=abs(Obs-Exp)
   
    yates=c(0.5,1,0.5)
    if (useYates==FALSE) yates=c(0,0,0)
    xtot=xtot-yates
    xtot=(xtot^2)/Exp
    xtot=colSums(xtot)
    pval=1-pchisq(xtot,1) # get p-value: high chi-square, low pvals
    
    # set chi-square and p-values to NA if INF or NaN - mostly homozygous
    xtot[xtot==Inf]=NA
    xtot[is.nan(xtot)]=NA
    pval[is.na(xtot)]=NA 
    
    # problem here - some values NA!!!!????
    HW=matrix(c(xtot,pval),curNumSNP,2)
    HW=HW[idindex,] 
    write(t(cbind(snps[firstSNP:lastSNP],HW)),paste(outputname,"hwstats.txt",sep="/"),ncolumns=3,sep="\t",append=T)  # write output
    
    firstSNP=lastSNP+1
    lastSNP=lastSNP+numHits 
  } 
  # end SNP loop
  
  sfStop() #turn off cluster
  
  rm(firstSNP,lastSNP,
  numIter,lastOnes,getSNP,curNumSNP,dbSNPQuery,
  hold,idindex,geno,gcscore,gcstats, misgeno,genotot,AAcount,ABcount,BBcount,alleletot,Bcount,
    Acount,q,p,snpcounts,snpfreqs,Exp,Obs,xtot,yates,pval,HW)
  
  cat("snp filtering complete.\n")
  print(Sys.time()-begin) # show runtime
}
# end snp QC



#################################
# slide descriptive statistics #
###############################
.sampleStats=function()
{
  begin=Sys.time() # keep track of how long it takes
  cat("running sample filtering...\n") 
  
  # slide summary
  sumslides=matrix(NA,numanim,4)
  rownames(sumslides)=animids
  colnames(sumslides)=c("Missing","AA","AB","BB")
  
  counter=0 
  
  numgeno=matrix(9,numsnp,numanim) # hold reshaped (numeric data)
  for (i in 1:numanim)
  {
    counter=counter+1
    if(counter==100) {counter=0;cat(paste(i,"\n"))}
    
    hold=dbGetQuery(dbcon,paste("select ",paste(fID,fsnp,fall1,fall2,fgc,sep=",")," from ",tSNP," where ",fID,"='",animids[i],"'",sep=""))
    names(hold)=c("ID","snp","allele1","allele2","gcscore")
    if (hasGC==FALSE) hold$gcscore=1 # set gcscore to one for all   
    
    geno=paste(as.character(hold$allele1),as.character(hold$allele2),sep="")
    geno[geno=="AA"]=0
    geno[geno=="AB"]=1
    geno[geno=="BA"]=1
    geno[geno=="BB"]=2
    geno[geno==paste(missingSymbol,missingSymbol,sep="")]=9
    
    geno[hold$gcscore<GCcutoff]=9 # change to 9 genotypes under GC score cutoff
    
    sumslides[i,1]=length(which(geno==9))
    sumslides[i,2]=length(which(geno==0))
    sumslides[i,3]=length(which(geno==1))
    sumslides[i,4]=length(which(geno==2))     
 
    geno=as.numeric(geno[match(snps,hold$snp)]) # order snps as per map position
    numgeno[,i]=geno   
  }      
  rownames(numgeno)=snps
  colnames(numgeno)=animids
  write.table(sumslides,paste(outputname,"samplecounts.txt",sep="/"),quote=F,sep="\t")
      
  cat("sample filtering complete.\n") 
  print(Sys.time()-begin) # show runtime
      
  return (numgeno)
}
# end sample QC


#################################
# sample correlations ##########
###############################
.sampleCorrelation=function()
{
  begin=Sys.time() # keep track of how long it takes  
  cat("running correlations...\n")
  
  # correlations
  animcor=cor(numgeno) # think about 9!!!
  animcor[which(is.na(animcor)==T)]=0 # if NA set to zero - will happen if whole genotype is the same e.g. all failed and got set to 9
  write.table(animcor,paste(outputname,"samplecorrelations.txt",sep="/"),quote=F,sep="\t") 
  cat("correlations complete.\n") 
  print(Sys.time()-begin) # show runtime  
}

