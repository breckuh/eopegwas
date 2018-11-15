.summaryAndPlots=function()
{
  # summary statistics - plots and tables
  # to stand alone only needs snparams
      
  require("gplots")
  require("gtools") 
  require("lattice")
  require("compiler")
  
  enableJIT(3) 
  
  begin=Sys.time() # keep track of how long it takes
  cat("building summaries and plots for report...\n")
  
  # read in data
  # map
  map=read.table(paste(outputname,"orderedmap.txt",sep="/"),header=T,sep="\t",colClasses=c("character","factor","numeric"),check.names=F)
  numsnp=nrow(map)
  
  # snp
  gcstats=read.table(paste(outputname,"gcstats.txt",sep="/"),header=T,sep="\t",colClasses=c("character",rep("numeric",12)),check.names=F)
  snpcounts=read.table(paste(outputname,"snpcounts.txt",sep="/"),header=T,sep="\t",colClasses=c("character",rep("numeric",7)))
  snpfreqs=read.table(paste(outputname,"snpfreqs.txt",sep="/"),header=T,sep="\t",colClasses=c("character",rep("numeric",7)))
  
  HW=read.table(paste(outputname,"hwstats.txt",sep="/"),header=T,sep="\t",colClasses=c("character",rep("numeric",2)))
  
  # samples
  sumslides = read.table(paste(outputname,"samplecounts.txt",sep="/"),header=T,sep="\t",colClasses=c("character",rep("numeric",4)),check.names=F)
  numanim=nrow(sumslides)
  
  animcor = read.table(paste(outputname,"samplecorrelations.txt",sep="/"),header=T,sep="\t",colClasses=c("character",rep("numeric",numanim)),check.names=F)
  animcor=as.matrix(animcor)  
  animcor[which(is.nan(animcor)==T)]=0 # set to zero missing correlations (e.g. samples with 9s across all SNPs)
  
  ######################### mapping summary and plots ####################################################
  # chromosome summary
  chromsum=summary(map$chrom)
  map$loc=as.numeric(map$loc)
  mapsum=NULL
  for (i in 1:length(names(chromsum)))
  {
    hold=map$loc[which(as.character(map$chrom)==names(chromsum)[i])]
    distance=sort(hold)[2:length(hold)]-sort(hold)[1:length(hold)-1]
    mapsum=rbind(mapsum,c(names(chromsum)[i],length(hold),(summary(distance))))
  }
  mapsum=mapsum[,-c(4,7)] # remove quartiles
  mapsum=data.frame(mapsum)
  names(mapsum)=c("chrom","num","min","median","mean","max")
  mapsum=mapsum[,c(1,2,3,6,5,4)]
  mapsum=mapsum[mixedorder(mapsum[,1]),] # sort chromosomes
  write.table(mapsum,paste(outputname,"mapsummary.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  # HW x mapping plot
  hold=data.frame(HW,map[,-1])
  pdf(file=paste(reportname,"HWmap.pdf",sep="/"),width=10,height=10) # save plot
  print(xyplot(chisquare~loc | chrom,data=hold,pch=".",col="black",xlab="location",ylab="HW chisquare",scales=list(x=list(alternating=4))))
  dev.off()
  
  # Median GC score x mapping plot
  hold=data.frame(gcstats,map[,-1])
  pdf(file=paste(reportname,"GCmap.pdf",sep="/"),width=10,height=10) # save plot
  print(xyplot(Median~loc | chrom,data=hold,pch=".",col="black",xlab="location",ylab="median GC score",scales=list(x=list(alternating=4))))
  dev.off()
  
  # MAF score x mapping plot
  # minor allele frequency
  # maf=apply(data.frame(A=snpcounts$A,B=snpcounts$B),1,function(x) if(is.na(x[1]) | x[1]<x[2]) x[1] else x[2])
  maf=apply(data.frame(A=snpcounts$A,B=snpcounts$B),1,function(x) if(x[1]<x[2]) x[1] else x[2])
  maf=maf/(snpcounts$A+snpcounts$B)
  hold=data.frame(snp=row.names(snpcounts),maf=maf,map[,-1])
  pdf(file=paste(reportname,"mafmap.pdf",sep="/"),width=10,height=10) # save plot
  print(xyplot(maf~loc | chrom,data=hold,pch=".",col="black",xlab="location",ylab="minor allele frequency",scales=list(x=list(alternating=4))))
  dev.off()    
  write.table(hold[,1:2],paste(outputname,"maf.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  # heterozygosity x mapping plot
  heterozygosity=snpcounts$AB/(snpcounts$AA+snpcounts$AB+snpcounts$BB) # heterozygosity
  hold=data.frame(snp=row.names(snpcounts),heterozygosity=heterozygosity,map[,-1])
  pdf(file=paste(reportname,"heteromap.pdf",sep="/"),width=10,height=10) # save plot
  print(xyplot(heterozygosity~loc | chrom,data=hold,pch=".",col="black",xlab="location",ylab="heterozygosity",scales=list(x=list(alternating=4))))
  dev.off()
  
  
  ############################ SNP summaries and Plots ###################################
  # call rates
  callrates=(snpcounts$AA+snpcounts$AB+snpcounts$BB)/numanim
  callsum=numeric(5)
  callsum[1]=length(which(callrates<0.9)) # number of calls under 0.9
  callsum[2]=length(which(callrates<0.95 & callrates>=0.9)) # number of calls under 0.95
  callsum[3]=length(which(callrates<0.99 & callrates>=0.95))# number of calls under 0.99
  callsum[4]=length(which(callrates<0.995 & callrates>=0.99)) # number of calls under 0.995
  callsum[5]=length(which(callrates>=0.995)) # number of calls under 1
  callsum=data.frame(rate=c("<0.9","0.9-0.95","0.95-0.99","0.99-0.995",">=0.995"),count=callsum,frequency=callsum/numsnp)
  write.table(callsum,paste(outputname,"callrates.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
   
  # call rates plot
  pdf(file=paste(reportname,"callrates.pdf",sep="/"),width=8,height=7)
  Yval=sort(callrates)
  forplot=length(Yval)
  Xval=1:length(Yval)/forplot*100
  plot(Xval,Yval,col="blue",xlab="proportion of snps",ylab="call rate frequency",main="distribution of call rates",
  cex.main=0.9,cex.lab=0.8,cex.axis=0.8,pch=".")
  abline(v=length(which(callrates<0.9))/forplot*100,col="black")
  abline(v=length(which(callrates<0.95))/forplot*100,col="red")
  abline(v=length(which(callrates<0.99))/forplot*100,col="green")
  abline(v=length(which(callrates<0.995))/forplot*100,col="orange")
  legend("bottomright",c("0.9","0.95","0.99","0.995"),cex=0.7,fill=c("black","red","green","orange"),bty="n")
  dev.off()
  
  # minor allele frequency
  forplot=seq(0,0.5,by=0.05)
  forplot2=numeric()
  for (i in 2:length(forplot)) forplot2[i-1]=length(which(maf<forplot[i] & maf>=forplot[i-1]))
  forplot2[length(forplot2)]=forplot2[length(forplot2)]+length(which(maf==0.5)) # add the 0.5s which were not included in loop
  
  pdf(file=paste(reportname,"maf.pdf",sep="/"),width=7,height=6)
  plot(forplot[-1]-0.025,forplot2,cex.main=0.9,cex.lab=0.8,cex.axis=0.8,xlab="minor allele frequency",ylab="number of snps",xlim=c(0,0.5),pch=20,col="blue",
   sub=paste("NAs: ",length(which(is.na(maf)==T)),"     MAF=0:", length(which(maf==0)),"     MAF<",snparams[7],": ",length(which(maf<snparams[7])),sep=""),cex.sub=0.8)
  lines(forplot[-1]-0.025,forplot2,col="blue")
  dev.off()
  
  # gc summary
  gssum=numeric()
  gssum[1]=sum(gcstats[,7])
  gssum[2]=sum(gcstats[,8])
  gssum[3]=sum(gcstats[,9])
  gssum[4]=sum(gcstats[,10])
  gssum[5]=sum(gcstats[,11])
  gssum=data.frame(band=names(gcstats)[7:11],count=gssum,frequency=gssum/(numsnp*numanim))
  write.table(gssum,paste(outputname,"gcbands.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  pdf(file=paste(reportname,"pieGC.pdf",sep="/"),width=8,height=8)
  pie(gssum[,2],col=1:5,labels=gssum[,1],main="distribution of GC scores",cex.main=0.9,cex.lab=0.8,cex.axis=0.8,
  xlab=paste("Missing:",sum(gcstats$Missing)))
  legend("topleft",paste(gssum[,1],": ",gssum[,2]," (",round(gssum[,3]*100,2),"%)",sep=""),fill=1:5,cex=0.7,bty="n")
  dev.off()
  
  gczero=length(which(gcstats$Max==0)) # number of snps with gc score = 0 in all chips
  gcnine=gcstats[,11]/numanim
  gcnine=length(which(gcnine>=0.9)) # number of snps where at least 90% of the samples have a gc score >0.9
  gccustom=1-gcstats[,12]/numanim
  gccustom=length(which(gccustom<(snparams[5]/100))) # number of snps where less than X% of the samples have a gc score > Y (user parameters)
  
  pdf(file=paste(reportname,"densityGC.pdf",sep="/"),width=8,height=8)
  plot(density(gcstats$Mean,na.rm=T),cex.main=0.9,cex.lab=0.8,cex.axis=0.8,main="mean density distribution of GC scores",col="blue")
  abline(v=GCcutoff,col="red")
  abline(v=0.9,col="red")
  legend("topleft",paste(c("100% GC=0:",paste(">",snparams[5],"% GC<",snparams[4],":",sep=""),">90% GC>0.9:"),c(gczero,gccustom,gcnine)),cex=0.7,bty="n")
  dev.off()
  
  # Hardy-Weinberg
  pdf(file=paste(reportname,"HW.pdf",sep="/"),width=10,height=8)
  Yval=sort(HW$pvalue)
  Xval=1:length(Yval)/length(Yval)*100
  plot(Xval,Yval,col="blue",pch=".",cex.main=0.9,cex.axis=0.8,cex.lab=0.8,cex.sub=0.8,xlab="proportion of SNPs", ylab="p-values",main="P-values for Hardy-Weinberg equilibrium",
  sub=paste("NAs: ",length(which(is.na(HW$pvalue)==T)),"     pval=0: ",length(which(HW$pvalue==0)),
  "     pval<",snparams[9],":",length(which(HW$pvalue<snparams[9])),
  sep=""))
  abline(v=length(which(HW$pvalue==0))/length(Yval)*100,col="black")
  abline(v=length(which(HW$pvalue<snparams[9]))/length(Yval)*100,col="red")
  legend("bottomright",c("pval=0",paste("pval<",snparams[9],sep="")),fill=c("black","red"),cex=0.8,bty="n")
  dev.off()
  
  # heterozygosity (He) and gene diversity (Ho)  
  genediv=1-((snpcounts$A/(snpcounts$A+snpcounts$B))^2+(snpcounts$B/(snpcounts$A+snpcounts$B))^2) # gene diversity
  
  up=mean(heterozygosity,na.rm=T)+snparams[8]*sd(heterozygosity,na.rm=T) # outliers 3 (user defined) SD
  down=mean(heterozygosity,na.rm=T)-snparams[8]*sd(heterozygosity,na.rm=T)
  hout=length(which(heterozygosity>up))
  hout=hout+length(which(heterozygosity<down))  # number of outliers
  holdSNPhetero=c(which(heterozygosity>up),which(heterozygosity<down))
  
  ymax=max(density(heterozygosity,na.rm=T)$y)
  ymax=max(c(ymax,density(genediv,na.rm=T)$y))
  
  pdf(file=paste(reportname,"heterozygosity.pdf",sep="/"),width=9,height=8)
  plot(NULL,xlim=c(0,1),ylim=c(0,ymax),col="blue",cex.main=0.9,cex.lab=0.8,cex.axis=0.8,xlab="density",ylab="frequency",
  main=paste("Heterozygosity (Ho) and gene diversity (He) density plot\nHo - mean:",round(mean(heterozygosity,na.rm=T),3),"     sd:",round(sd(heterozygosity,na.rm=T),3),
  "\nHe - mean:",round(mean(genediv,na.rm=T),3),"     sd:",round(sd(genediv,na.rm=T),3)),
  sub=paste("mean: black line  /  ",snparams[8],"SD: red line  /  number of outliers: ",hout,sep=""),cex.sub=0.8)
  
  lines(density(heterozygosity,na.rm=T),lty=1,col="blue")
  lines(density(genediv,na.rm=T),lty=2,col="green")
  
  abline(v=mean(heterozygosity,na.rm=T))
  abline(v=up,col="red")
  abline(v=down,col="red")
  legend("topright",c("Ho","He"),lty=c(1,2),col=c("blue","green"),bty="n",cex=0.8)
  dev.off()
  
  ######################### Sample plots (summaries generated in preprocessing) #############################
  # heatmap of correlations and correlation plot
  hmcol=greenred(256)
  pdf(file=paste(reportname,"heatcorrelation.pdf",sep="/"),width=10,height=10)
  heatmap(animcor,col=hmcol,symm=T,labRow=" ",labCol=" ")
  dev.off()
  
  #dev.print(file=paste(reportname,"heatcorrelation.pdf",sep="/"),device=pdf,width=10,height=10) # save plot
  png(file=paste(reportname,"heatcorrelation.png",sep="/"),width=1024,height=768)
  heatmap(animcor,col=hmcol,symm=T,labRow=" ",labCol=" ")
  dev.off()
  
  # get similar animals
  updiag=upper.tri(animcor,diag=F)
  hold=updiag
  hold[updiag]=animcor[updiag]
  rownames(hold)=rownames(animcor)
  colnames(hold)=colnames(animcor)
  similar=matrix(NA,1,3)
  j=1
  for (i in 1: length(hold[,1]))
  {
    index=which(hold[,i]>sampleparams[3])
    if (length(index)>0)
    {
      temp=matrix(NA,length(index),3)
      for (k in 1:length(index))
      {
        temp[k,1]=rownames(hold)[index[k]]
        temp[k,2]=colnames(hold)[i]
        temp[k,3]=round(hold[index[k],i],3)
      }
      similar=rbind(similar,temp)
      j=j+1
    }
  }
  similar=similar[-1,]
  similar=similar[order(similar[,1],similar[,2]),]
  similar=data.frame(similar)
  names(similar)=c("sample1","sample2","correlation")
  write.table(similar,paste(outputname,"similar.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  pdf(file=paste(reportname,"correlation.pdf",sep="/"),width=10,height=8) # save plot
  plot(density(animcor[updiag]),col="blue",cex.main=0.9,cex.axis=0.8,cex.lab=0.8,ylab="correlation",xlab="",main="Correlation between samples",pch=".")
  legend("topleft",cex=0.8,bty="n",
  paste(c("min: ", "max: ", "mean: ", "median: ",">0.9: ","<0.1: "), c(round(min(animcor[updiag]),3),round(max(animcor[updiag]),3),round(mean(animcor[updiag]),3),
  round(median(animcor[updiag]),3),length(which(animcor[updiag]>0.9)),length(which(animcor[updiag]< 0.1)))))
  dev.off()
  
  # sample summary
  hold=rowSums(sumslides[,2:4])
  hold=data.frame(sample=names(hold),calls=hold,frequency=hold/numsnp)
  samplecalls=hold
  sumslides2=numeric()
  sumslides2[1]=min(hold$frequency,na.rm=T)
  sumslides2[2]=max(hold$frequency,na.rm=T)
  sumslides2[3]=mean(hold$frequency,na.rm=T)
  sumslides2[4]=length(which(hold$frequency<0.97))
  sumslides2[5]=length(which(hold$frequency<sampleparams[1]))
  sumslides2=data.frame(statistic=c("min","max","mean","<0.97",paste("<",sampleparams[1],sep="")),value=sumslides2)
  write.table(sumslides2,paste(outputname,"samplestats.txt",sep="/"),quote=F,sep="\t",row.names=F)
  holdanimcalls=row.names(sumslides)[which(hold$frequency<sampleparams[1])]
  
 
  # slide heterozygosity
  animids=rownames(sumslides)
  samplehetero=sumslides[,3]/(sumslides[,2]+sumslides[,3]+sumslides[,4])
  up=mean(samplehetero,na.rm=T)+sampleparams[2]*sd(samplehetero,na.rm=T) # outliers 3 SD
  down=mean(samplehetero,na.rm=T)-sampleparams[2]*sd(samplehetero,na.rm=T)
  hsout=length(which(samplehetero>up))
  hsout=hsout+length(which(samplehetero<down))  # number of outliers
  holdanimhetero=c(which(samplehetero>up),which(samplehetero<down))
  holdanimhetero=animids[holdanimhetero]     
  
  sd1=mean(samplehetero,na.rm=T)-sampleparams[2]*sd(samplehetero,na.rm=T)
  sd2=mean(samplehetero,na.rm=T)+sampleparams[2]*sd(samplehetero,na.rm=T)
  xlims=c(min(samplehetero,na.rm=T),max(samplehetero,na.rm=T))
  if (xlims[1]>sd1) xlims[1]=sd1-sd1*0.01 else xlims[1]=xlims[1]+xlims[1]*0.01
  if (xlims[2]<sd2) xlims[2]=sd2+sd2*0.01 else xlims[2]=xlims[2]+xlims[2]*0.01
  pdf(file=paste(reportname,"samplehetero.pdf",sep="/"),width=6,height=6) # save plot
  hold=sort(samplehetero)
  plot(hold,1:length(hold),col="blue",cex.main=0.9,cex.axis=0.8,cex.lab=0.8,ylab="sample",xlab="heterozygosity", xlim=xlims,
  main=paste("Sample heterozygosity\nmean:",round(mean(samplehetero,na.rm=T),3),"    sd:",round(sd(samplehetero,na.rm=T),3)),
  sub=paste("mean: black line     ",sampleparams[2],"SD: red line     number of outliers:",hsout),cex.sub=0.8)
  abline(v=mean(samplehetero,na.rm=T))
  abline(v=mean(samplehetero,na.rm=T)-sampleparams[2]*sd(samplehetero,na.rm=T),col="red")
  abline(v=mean(samplehetero,na.rm=T)+sampleparams[2]*sd(samplehetero,na.rm=T),col="red")
  dev.off()
  
   
  
  ######################## good/bad plots ######################################## 
  dbcon<<-dbConnect(dbDriver("SQLite"), dbname = dbname)
  source("GetTablesFields.r")  
  
  for (i in 1:2)
  {
    if (i==1)  index=which(gcstats$Mean==min(gcstats$Mean,na.rm=T))[1] # bad slide index
    if (i==2)  index=which(gcstats$Mean==max(gcstats$Mean,na.rm=T))[1] # good slide index
    
    getSNP=row.names(gcstats)[index]
    
    hold=dbGetQuery(dbcon,paste("select * from ", tSNP," where ",fsnp,"='", rownames(gcstats)[index],"'",sep=""))
    if (hasXY==TRUE)  
    {
      hold=hold[,c(fID,fsnp,fall1,fall2,fgc,fx,fy)]
      names(hold)=c("ID","snp","allele1","allele2","gcscore","x","y") 
    } else  
    {
      hold=hold[,c(fID,fsnp,fall1,fall2,fgc)]
      names(hold)=c("ID","snp","allele1","allele2","gcscore")
    }
    if (hasGC==FALSE) hold$gcscore=1 
  
    tempgc=as.numeric(as.character(hold$gcscore))
    tempgc[is.na(tempgc)==T]=0 # set NAs to zero
    tempgc[is.nan(tempgc)==T]=0 # set NANs to zero
    
   geno=paste(hold$allele1,hold$allele2,sep="")
      geno[geno=="AA"]=0
      geno[geno=="AB"]=1
      geno[geno=="BA"]=1
      geno[geno=="BB"]=2
      geno[geno==paste(missingSymbol,missingSymbol,sep="")]=NA
    
    # setup screen
    if (i==1) pdf(file=paste(reportname,"bad.pdf",sep="/"),width=10,height=8) # just to save a bad example for the report
    if (i==2) pdf(file=paste(reportname,"good.pdf",sep="/"),width=10,height=8) # just to save a good example for the report
    
    par(mfrow = c(2, 2), oma = c(2, 2, 2, 2))  
    # calls plot
    if(hasXY==TRUE)
    {
      geno=as.factor(geno)
      plot(hold$x,hold$y,col=as.numeric(geno),pch=as.numeric(geno),xlab="x",ylab="y",
      main=hold$snp[1],cex.main=0.9,cex.lab=0.9,cex.axis=0.8) 
            
      legend("bottomleft",paste(" (",summary(geno),")",sep=""),
      col= 1:length(levels(geno)),
      pch= 1:length(levels(geno)),cex=0.7)    
    } else
    {
     plot(NA,NA,cex.main=0.9,cex.lab=0.9,cex.axis=0.8,xlim=c(0,1),ylim=c(0,1),
     xlab="x",ylab="y", main=hold$snp[1])   
    }
  
    # gc score plot    
    plot(tempgc, ylim=c(0,1),xlab="samples",ylab="gc score",cex.lab=0.8,cex.axis=0.8,col="blue",
    sub=paste("NA or 0:",gcstats[index,1],"     gc<",GCcutoff,":",gcstats[index,12],"     gc>=",GCcutoff,":",(numanim-gcstats[index,12]),sep=""),cex.sub=0.8)
    abline(h=GCcutoff)
  
    # alleles   
    Obs=as.numeric(snpfreqs[index,1:3])
    names(Obs)=c("MIS","A","B")
    barplot(as.numeric(snpfreqs[index,1:3]),col=c("black","red","blue"),cex.axis=0.8,cex.lab=0.8,cex.sub=0.8,cex.main=0.8,main="allelic frequencies",
    sub=paste("-: ",round(snpfreqs[index,1],2),"     A: ",round(snpfreqs[index,2],2),"     B: ",round(snpfreqs[index,3],2),sep=""))
  
    # genotypes plot   
    Obs=as.numeric(snpcounts[index,5:7])
    names(Obs)=c("AA","AB","BB")
  
    p=snpfreqs[index,2]
    q=snpfreqs[index,3]
    Exp=c(p^2,2*p*q,q^2)
    Exp=Exp*sum(snpcounts[index,5:7])
    names(Exp)=c("AA","AB","BB")     
      
    barplot(c(MIS=snpcounts[index,4],Exp,Obs),col=c("black","red2","lightgreen","lightblue","red","green","blue"),cex.axis=0.8,cex.lab=0.8,cex.sub=0.8,cex.main=0.8,
    main=paste("genotypic frequencies\nHW p-value:",signif(HW$pvalue[index],4)),
    sub=paste("MIS=missing, left expected, right observed"))
    abline(h=0)
       
    dev.off()  
  }
  rm(hold,Exp,Obs,tempgc,p,q,i,geno,getSNP) # clean up
  dbDisconnect(dbcon) # close DB connection
    
  cat("summaries and plots complete.\n")
  print(Sys.time()-begin) # show runtime
  
  #####################################
  ## Get SNPs and samples to remove ##
  ###################################
  begin=Sys.time() # keep track of how long it takes
  cat("building output files...\n")
  
  # summary of numbers removed/filtered
  QCcrit=c(paste(">",snparams[1]," percent genotyping fail",sep=""),
  paste("median GC scores <",snparams[2],sep=""),
  paste("all GC scores",snparams[3]),
  paste("GC <",snparams[4]," in less than ",snparams[5]," percent samples",sep=""),
  paste("100 percent homozygous"),
  paste("MAF <",snparams[7],sep=""),
  paste("heterozygosity ",snparams[8],"SD",sep=""),
  paste("Hardy-Weinberg at",snparams[9])
  )
  
  #### maf reports numbers of homozygous SNP but does not necessarily remove unless snparams[6] set to zero.
  QCnum=c(
  length(which(gcstats$Missing>(numanim/100*snparams[1]))),
  length(which(gcstats$Median<snparams[2])),
  gczero,
  gccustom,
  length(which(maf==0)),
  length(which(maf>0 & maf<snparams[7])),
  hout,
  length(which(HW$pvalue<snparams[9]))
  )
  reject=data.frame("SNP criteria"=QCcrit,number=QCnum,check.names=F)
  write.table(reject,paste(outputname,"rejectSNPsummary.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  QCcrit=c(paste("call rates <",sampleparams[1],sep=""),
  paste("correlation >",sampleparams[3],sep=""),
  paste("heterozygosity ",sampleparams[2],"SD",sep="")
  )
  QCnum=c(sumslides2$value[5],
  length(which(animcor[updiag]>sampleparams[3])),
  hsout
  )
  reject2=data.frame("sample criteria"=QCcrit,number=QCnum,check.names=F)
  write.table(reject2,paste(outputname,"rejectSamplesummary.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  mapexcludeSNP=intersect(names(chromsum),mapexcludeSNP) # check if excluded snps actually exist in snp map
  QCnum=numeric()
  QCcrit=character()
  unmapped=0
  if (length(mapexcludeSNP)>0)
  {
    for (i in 1:length(mapexcludeSNP))
    {
      QCcrit=c(QCcrit,paste("Chromosome",mapsum$chrom[which(mapsum$chrom==mapexcludeSNP[i])]))
      unhold=as.numeric(as.character(mapsum$num[which(mapsum$chrom==mapexcludeSNP[i])]))
      QCnum=c(QCnum,unhold)
      unmapped=unmapped+unhold  
    }
  } else
  {
    QCcrit="No chromosomes excluded"
    QCnum=0
  }
  rm(unhold)
  reject3=data.frame("mapping criteria"=QCcrit,number=QCnum,check.names=F)
  write.table(reject3,paste(outputname,"rejectMapsummary.txt",sep="/"),quote=F,sep="\t",row.names=F)
  
  # indexes/names of SNPs to remove
  snpindex=which(gcstats$Missing>(numanim/100*snparams[1]))
  snpindex=c(snpindex,which(gcstats$Median<snparams[2]))
  snpindex=c(snpindex,which(gcstats$Max==0))
  
  hold=1-gcstats[,12]/numanim
  snpindex=c(snpindex,which(hold<(snparams[5]/100))) # number of snps where less than X% of the samples have a gc score > Y (user parameters)
  
  # maf
  snpindex=c(snpindex,which(maf<=snparams[6]))
  snpindex=c(snpindex,which(maf>0 & maf<snparams[7]))
  
  # heterozygosity
  snpindex=c(snpindex,holdSNPhetero)
  
  # HW
  snpindex=c(snpindex,which(HW$pvalue<snparams[9]))
  
  # join into unique list
  snpindex=c(row.names(gcstats)[snpindex])
  snpindex=unique(snpindex)
  write.table(snpindex,paste(outputname,"rejectSNPnames.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
  
  # indexes/names of slides to remove
  animindex=c(holdanimcalls,holdanimhetero) # calls and heterozygosity
  animindex=unique(animindex) # unique list
  write.table(animindex,paste(outputname,"rejectSamplenames.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
  
  # indexes/names of snps to remove due to mapping - user defined
  map$snp=as.character(map$snp)
  mapindex=character()
  if (length(mapexcludeSNP)>0)
  {
    for (i in 1:length(mapexcludeSNP))
    {
      mapindex=c(mapindex,map$snp[which(map$chrom==mapexcludeSNP[i])])
    }
    mapindex=unique(mapindex)
  }
  write.table(mapindex,paste(outputname,"rejectMapnames.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
  write.table(unique(c(mapindex,snpindex)),paste(outputname,"rejectMapAndSNPs.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
  
  # build and save final reshaped data - only GC rejected set to 9 or snp/samples rejected set to 9 (or excluded)
  # change symbols if defined by user
  if (exists("numgeno")==TRUE) # check if numgeno exists in run
  {    
    numgeno2=numgeno
    if (length(animindex>0)) numgeno2[,match(animindex,colnames(numgeno))]=9
    numgeno2[match(unique(c(mapindex,snpindex)),rownames(numgeno)),]=9
    
    numexc=length(which(numgeno2==9)) # total number of snps excluded (gc + whole snp/sample)
    write.table(numexc,paste(outputname,"numExcluded.txt",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
    
    if (missingReshaped!=9)id1=which(numgeno==9)
    if (genoReshaped[1]!=0)id2=which(numgeno==0)
    if (genoReshaped[2]!=1)id3=which(numgeno==1)
    if (genoReshaped[3]!=2)id4=which(numgeno==2)
    
    if (excludeReshaped==TRUE) # remove all rejected SNP and samples
    {
      if (missingReshaped!=9) numgeno[id1]=missingReshaped
      if (genoReshaped[1]!=0) numgeno[id2]=genoReshaped[1]
      if (genoReshaped[2]!=1) numgeno[id3]=genoReshaped[2]
      if (genoReshaped[3]!=2) numgeno[id4]=genoReshaped[3]
    
      idx1=match(animindex,colnames(numgeno))
      idx2=match(unique(c(mapindex,snpindex)),rownames(numgeno))
      if (length(idx1>0)) numgeno=numgeno[,-idx1]
      if (length(idx2>0))numgeno=numgeno[-idx2,]     
    
      write.table(data.frame(animalID=colnames(numgeno)),paste(outputname,"animalID.txt",sep="/"),quote=F,row.names=F,sep="\t") # save animal id order in case needed to use flat files
      write.table(data.frame(SNP=rownames(numgeno)),paste(outputname,"snpID.txt",sep="/"),quote=F,row.names=F,sep="\t") # save snp id order in case needed to use flat files
      write.table(numgeno,paste(outputname,"reshaped.txt",sep="/"),row.names=F,col.names=F,quote=F,sep=" ") # save file in format for downstream analysis
    } else if (excludeReshaped==FALSE & setReshaped9==TRUE) # keep all but set to missing value
    {
      if (missingReshaped!=9) numgeno2[id1]=missingReshaped
      if (genoReshaped[1]!=0) numgeno2[id2]=genoReshaped[1]
      if (genoReshaped[2]!=1) numgeno2[id3]=genoReshaped[2]
      if (genoReshaped[3]!=2) numgeno2[id4]=genoReshaped[3]
    
      write.table(data.frame(animalID=colnames(numgeno2)),paste(outputname,"animalID.txt",sep="/"),quote=F,row.names=F,sep="\t") # save animal id order in case needed to use flat files
      write.table(data.frame(SNP=rownames(numgeno2)),paste(outputname,"snpID.txt",sep="/"),quote=F,row.names=F,sep="\t") # save snp id order in case needed to use flat files
      write.table(numgeno2,paste(outputname,"reshaped.txt",sep="/"),row.names=F,col.names=F,quote=F,sep=" ") # save file in format for downstream analysis
    } else
    {
      if (missingReshaped!=9) numgeno[id1]=missingReshaped
      if (genoReshaped[1]!=0) numgeno[id2]=genoReshaped[1]
      if (genoReshaped[2]!=1) numgeno[id3]=genoReshaped[2]
      if (genoReshaped[3]!=2) numgeno[id4]=genoReshaped[3]
    
      write.table(data.frame(animalID=colnames(numgeno)),paste(outputname,"animalID.txt",sep="/"),quote=F,row.names=F,sep="\t") # save animal id order in case needed to use flat files
      write.table(data.frame(SNP=rownames(numgeno)),paste(outputname,"snpID.txt",sep="/"),quote=F,row.names=F,sep="\t") # save snp id order in case needed to use flat files
      write.table(numgeno,paste(outputname,"reshaped.txt",sep="/"),row.names=F,col.names=F,quote=F,sep=" ") # save file in format for downstream analysis
    }
  } 
  else numexc=readLines(paste(outputname,"numExcluded.txt",sep="/")) # just for report
  
  ########################################################
  ## prepare flat files for adding results to database ###
  #######################################################
  
  # build snpfilter file for database - results of snp QC
  snpfilter=data.frame(snp=rownames(snpcounts),good=rep("Y",length(maf)),maf=maf,heterozygosity=heterozygosity,HW,callrates,snpcounts,mapped=rep("Y",length(maf)))
  snpfilter$mapped=as.character(snpfilter$mapped)
  snpfilter$good=as.character(snpfilter$good)
  snpfilter$mapped[match(mapindex,snpfilter$snp)]="N"  # set map excluded to N
  snpfilter$good[match(snpindex,snpfilter$snp)]="N" # set rejected snp to N
  write.table(snpfilter,paste(outputname,"snpfilter.txt",sep="/"),quote=F,row.names=F,col.names=T,sep="\t")
  
  # build samplefilter file for datbase - results of sample QC
  samplefilter=data.frame(id=animids,good=rep("Y",length(samplehetero)),heterozygosity=samplehetero,callrates=samplecalls$frequency,sumslides,check.names=F)
  samplefilter$good=as.character(samplefilter$good)
  samplefilter$good[match(animindex,samplefilter$id)]="N" # set rejected samples to N
  write.table(samplefilter,paste(outputname,"samplefilter.txt",sep="/"),quote=F,row.names=F,col.names=T,sep="\t")
  
  ##################################################
  ######### write variable values for report ######
  ################################################
  
  Lnumgeno=numsnp*numanim
  
  repdata=c(numanim,
  numsnp,
  length(animindex),
  length(unique(c(mapindex,snpindex))),
  length(which(callrates<((100-snparams[1])/100))),
  snparams[1],
  GCcutoff,
  gczero,
  gcnine,
  gccustom,
  snparams[5],
  round(mean(gcstats$Mean,na.rm=T),3),
  round(mean(gcstats$Median,na.rm=T),3),
  length(which(maf==0)),
  length(which(maf>0 & maf<snparams[7])),
  snparams[7],
  round(mean(heterozygosity,na.rm=T),3),
  round(sd(heterozygosity,na.rm=T),3),
  hout,
  length(which(is.na(HW$pvalue)==T)),
  length(which(HW$pvalue==0)),
  length(which(HW$pvalue<snparams[9]/numsnp)),
  length(which(HW$pvalue<snparams[9])),
  snparams[9],
  sampleparams[1],
  round(mean(animcor[updiag],na.rm=T),3),
  round(min(animcor[updiag],na.rm=T),3),
  round(max(animcor[updiag],na.rm=T),3),
  length(which(animcor[updiag]>sampleparams[3])),
  sampleparams[3],
  round(mean(samplehetero,na.rm=T),3),
  round(sd(samplehetero,na.rm=T),3),
  sampleparams[2],
  hsout,
  unmapped,
  Lnumgeno,
  numexc
  )
  write.table(data.frame(repdata),paste(outputname,"ForReport.txt",sep="/"),quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  cat(paste("output files saved to directory",outputname,"\n"))
  print(Sys.time()-begin) # show runtime   
}






