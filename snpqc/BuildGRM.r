if (buildGRM==T)
{
  cat("reading files for GRM...\n")
  # snp and samples - always read in from file; not efficient but easy to get right
  snps=as.vector(read.table(paste(outputname,"snpID.txt",sep="/"),sep="\t",header=T)[,1])
  animids=as.vector(read.table(paste(outputname,"animalID.txt",sep="/"),sep="\t",header=T)[,1])

  rejectsnps=as.vector(read.table(paste(outputname,"rejectMapAndSNPs.txt",sep="/"),sep="\t",header=F)[,1])
  rejectsamples=readLines(paste(outputname,"rejectSamplenames.txt",sep="/"))

  # match rejected with sample and snp ids
  manim=match(rejectsamples,animids)
  manim=manim[which(is.na(manim)==F)] # check if not already excluded
  
  msnp=match(rejectsnps,snps)
  msnp=msnp[which(is.na(msnp)==F)]

  # check if numgeno exists in run
  if (exists("numgeno")==FALSE) numgeno=matrix(scan(paste(outputname,"reshaped.txt",sep="/"),sep=" "),length(snps),length(animids),byrow=T)

  if (length(msnp)>0) numgeno=numgeno[-msnp,] # remove excluded snps
  if (length(manim)>0)
  {
    numgeno=numgeno[,-manim] # remove excluded samples
    write.table(animids[-manim],paste(outputname,"GRMsampleIDs.txt",sep="/"),sep="\t",quote=F,col.names=F,row.names=F)
  }
  else write.table(animids,paste(outputname,"GRMsampleIDs.txt",sep="/"),sep="\t",quote=F,col.names=F,row.names=F)
  #rm(snps,animids,rejectsnps,rejectsamples,manim,msnp) # cleanup
  
  cat("filtering data...\n")
  # filter
  numnine=apply(numgeno,1,function (x) length(which(x==missingReshaped)))

  compSNP=which(numnine<mean(numnine))
  numgeno=numgeno[compSNP,] # should usually give a large number of snps

  fillin = function(x) # ok if population reasonably homogeneous and no bias in genotyping problems
  {
    avg=mean(x[which(x!=missingReshaped)])
    x[which(x==missingReshaped)]=avg
    return (x)
  }
  for (i in 1:nrow(numgeno)) numgeno[i,]=fillin(numgeno[i,]) # put average value of genotypes in missing genotypes - very slow but less memory intensive

  #rm(numnine,compSNP)

  cat("building GRM...\n")
  p=apply (numgeno,1,function(x) sum(x)/(length(x)*2)) # frequency of second (minor) allele = num alleles / total num alleles

  # if 2 is not minor allele can replace 0 with 2 and 2 with 0 (but will make no difference to GRM)
  #MAFrep=which(p>0.5)
  #if (length(MAFrep)>0)
  #{
  #  MAFswap=numgeno[MAFrep,]

  #  index1=which(MAFswap==0)
  #  index2=which(MAFswap==2)

  #  MAFswap[index1]=2
  #  MAFswap[index2]=0

  #  numgeno[MAFrep,]=MAFswap
  #  p=apply (numgeno,1,function(x) sum(x)/(length(x)*2)) # frequency of second (minor) allele = num alleles / total num alleles
  #}

  numgeno=numgeno-1 # recode matrix as -1, 0, 1

  # GOF - scales G to analogous of NRM
  P=2*(p-0.5) # deviation from 0.5 - P should use base population frequencies!
  Z=numgeno-P
  ZtZ = t(Z) %*% Z
  d=2*sum(p*(1-p))
  G=ZtZ/d
  
  #rm(P,p,Z,ZtZ,d,numgeno)

  write.table(G,paste(outputname,"GRM.txt",sep="/"),sep=" ",quote=F,col.names=F,row.names=F)

  cat("GRM complete.\n")

  # singular value decomposition of GRM
  #SVD=svd(G)
  #plot(SVD$v[,1],SVD$v[,2],cex.main=0.9,main="Singular Value Decomposition",pch=16,xlab="PC1",ylab="PC2")
  
  #heatmap(G,symm=T,col=gray.colors(16,start=0,end=1))
  #dev.print("SVD12.pdf",device=pdf,width=10,height=8)
}


