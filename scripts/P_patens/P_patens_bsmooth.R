if (!require(bsseq)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("bsseq")
}

library(bsseq)

fitData <- function(dfMethData,m_context){
  ## In order to read the data, the raw data have to be adjusted
  # chromosome name, strand and context is missing
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
  
  # Smooth methylation data
  dataBSmooth <- BSmooth(m_dataBS, mc.cores=1, verbose=T)
  
  pData(m_dataBS)$Pair <- c("replicate1","replicate2","replicate3")
  pData(m_dataBS)$Type <- c("test","test","test","control","control","control")
  
  dataBS_cov <- getCoverage(dataBSmooth ==10 )
  
  
  BS_cov_filtered <- which(rowSums(dataBS_cov[,m_dataBS$Type == "test"]>= 2)>=2 & 
                             rowSums(dataBS_cov[,m_dataBS$Type == "control"]>=2)>=2)
  
  dataBSmooth <- dataBSmooth[BS_cov_filtered,]
  
  dataBSmooth_tstat <- BSmooth.tstat(dataBSmooth, 
                                     group1 = c("percentC_sample1", "percentC_sample2","percentC_sample3"),
                                     group2 = c("percentC_sample4", "percentC_sample5","percentC_sample6"), 
                                     estimate.var="same", local.correct = T, verbose = T)
  
  # finding dmrs
  dmrs <- dmrFinder(dataBSmooth_tstat,maxGap=400)
  dmrs_sub <- subset(dmrs , n >= 3 & abs(meanDiff) >= 0.1)
  dmrs_sub_ordered <- dmrs_sub[order(dmrs_sub$start),]
  
  # extracting dmrs for bed file
  #dmrs_for_bed <- dmrs_sub[order(dmrs_sub[1:3]$start),]#,1:3] for the first 3 columns
  write.table(dmrs_sub_ordered[1:3],paste("output_dir",newFolder,"/bsmooth_dmrs_",newFolder,".bed",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
  write.table(dmrs_sub_ordered,paste("output_dir",newFolder,"/bsmooth_dmrs_with_info_",newFolder,".bed",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
}
##############################################################
setwd("input_dir")

sim_cpg_g129 <- read.delim("/data/G129_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_cpg_r143 <- read.delim("/data/R143_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_r143 <- read.delim("/data/R143_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_g129 <- read.delim("/data/G129_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_r143 <- read.delim("/data/R143_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))      
sim_chh_g129 <- read.delim("/data/G129_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))                  


dataBS_cpg_g129 <- fitData(sim_cpg_g129,"CpG")
dataBS_cpg_r143 <- fitData(sim_cpg_r143,"CpG")
dataBS_chg_g129 <- fitData(sim_chg_g129,"CHG")
dataBS_chg_r143 <- fitData(sim_chg_r143,"CHG")
dataBS_chh_r143 <- fitData(sim_chh_r143,"CHH")
dataBS_chh_g129 <- fitData(sim_chh_g129,"CHH")



peakRAM(DMR_calling(dataBS_cpg_g129,"G129", "CpG"))
peakRAM(DMR_calling(dataBS_cpg_r143,"R143", "CpG"))
peakRAM(DMR_calling(dataBS_chg_r143,"R143", "CHG"))
peakRAM(DMR_calling(dataBS_chg_r143,"G129", "CHG"))
peakRAM(DMR_calling(dataBS_chh_r143,"R143", "CHH"))
peakRAM(DMR_calling(dataBS_chh_g129,"G129", "CHH"))


