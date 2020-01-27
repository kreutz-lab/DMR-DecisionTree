if (!require(bsseq)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("bsseq")
}

library(bsseq)

fitData <- function(dfMethData,m_context){
  methDf <- cbind(c(rep("chr1",dim(dfMethData)[1])),
                  c(rep("+",dim(dfMethData)[1])),
                  c(rep(m_context,dim(dfMethData)[1])),
                  dfMethData,stringsAsFactors=F)
  colnames(methDf)[1:4] <- c("chr","strand","context","pos")
  
  #change the proportinal thytosine values to absolute values
  for (i in 2:(length(methDf)/2-1)){
    methDf[i*2+2] <- methDf[i*2+1]*((methDf[i*2+2]/100))
  }
  
  # Creates a bsseq-class from Bsmooth
  BS <- BSseq(pos = methDf$pos, chr = methDf$chr, 
              M  = as.matrix(methDf[seq(6,length(methDf),by=2)]),
              Cov = as.matrix(methDf[seq(5,length(methDf),by=2)]))
  return(BS)
}


DMR_calling <- function(m_dataBS,m_ecotype, m_context){
  newFolder <- paste(m_ecotype,"_",m_context,sep="")
  if (!dir.exists(paste("output_dir",newFolder,sep=""))){
    dir.create(paste("output_dir",newFolder,sep=""))
  }
  
  dataBSmooth <- BSmooth(m_dataBS, mc.cores=1, verbose=T)
  
  pData(m_dataBS)$Pair <- c("replicate1","replicate2","replicate3")
  pData(m_dataBS)$Type <- c("test","test","test","control","control","control")
  
  dataBS_cov <- getCoverage(dataBSmooth ==10 )
  
  BS_cov_filtered <- which(rowSums(dataBS_cov[,m_dataBS$Type == "test"]>= 2)>=2 & 
                             rowSums(dataBS_cov[,m_dataBS$Type == "control"]>=2)>=2)
  
  dataBSmooth <- dataBSmooth[BS_cov_filtered,]
  
  dataBSmooth_tstat <-  BSmooth.tstat(dataBSmooth, 
                                      group1 = c("percentC_sample1", "percentC_sample2","percentC_sample3"),
                                      group2 = c("percentC_sample4", "percentC_sample5","percentC_sample6"), 
                                      estimate.var="same", local.correct = T, verbose = T)
  
  write.table(dmrs_sub_ordered[1:3],paste("output_dir",newFolder,"/bsmooth_dmrs_",newFolder,".bed",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
  write.table(dmrs_sub_ordered,paste("output_dir",newFolder,"/bsmooth_dmrs_with_info_",newFolder,".bed",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
}

setwd("input_dir/")

sim_arabidopsis_chg <- read.delim("input_dir/arabidopsis_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_arabidopsis_chh <- read.delim("input_dir/arabidopsis_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_arabidopsis_cpg <- read.delim("input_dir/arabidopsis_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))


dataBS_chg_arabidopsis <- fitData(sim_arabidopsis_chg,"CHG")
dataBS_chh_arabidopsis <- fitData(sim_arabidopsis_chh,"CHH")
dataBS_cpg_arabidopsis <- fitData(sim_arabidopsis_cpg,"CpG")

peakRAM(DMR_calling(dataBS_chg_arabidopsis,"arabidopsis", "CHG"))
peakRAM(DMR_calling(dataBS_chh_arabidopsis,"arabidopsis", "CHH"))
peakRAM(DMR_calling(dataBS_cpg_arabidopsis,"arabidopsis", "CpG"))

