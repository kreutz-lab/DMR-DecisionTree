
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
  
  for (i in 2:(length(methDf)/2-1)){
    methDf[i*2+2] <- methDf[i*2+1]*((methDf[i*2+2]/100))
  }
  
  BS <- BSseq(pos = methDf$pos, chr = methDf$chr, 
              M  = as.matrix(methDf[seq(6,length(methDf),by=2)]),
              Cov = as.matrix(methDf[seq(5,length(methDf),by=2)]))
  return(BS)
}

DMR_calling <- function(m_dataBS,m_ecotype, m_context){
  newFolder <- paste(m_ecotype,"_",m_context,sep="")
  if (!dir.exists(paste("output/bsmooth/",newFolder,sep=""))){
    dir.create(paste("output/bsmooth/",newFolder,sep=""))
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
  
  # Plot of marginal distribution of the t-statistics corrected and uncorrected
  pdf(paste("output/bsmooth/",newFolder,"/distribution_t-stat.pdf",sep=""))
  plot(dataBSmooth_tstat)
  dev.off()
  
  
  # plotting the methylation by myself
  cpgs_meth <- matrix(dataBSmooth@rowRanges@ranges@start)
  
  for (i in 1:6){
    cpgs_meth <- cbind(cpgs_meth,round(as.vector(getMeth(dataBSmooth,type="smooth")[,i]),5))
  }
  
  options("scipen"=50)
  write.table(cpgs_meth,paste("output/bsmooth/",newFolder,"/methylation_levels_smooth_",newFolder,".txt",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
  
  # finding dmrs
  dmrs <- dmrFinder(dataBSmooth_tstat,maxGap=400)
  dmrs_sub <- subset(dmrs , n >= 3 & abs(meanDiff) >= 0.2)
  dmrs_sub_ordered <- dmrs_sub[order(dmrs_sub$start),]
  
  write.table(dmrs_sub_ordered[1:3],paste("output/bsmooth/",newFolder,"/bsmooth_dmrs_",newFolder,".bed",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
  write.table(dmrs_sub_ordered,paste("output/bsmooth/",newFolder,"/bsmooth_dmrs_with_info_",newFolder,".bed",sep=""),
              row.names=F,col.names=F,sep="\t",quote=F)
  
}

pData <- pData(dataBSmooth)
pData$col <- rep(c("red", "blue"),each = 3)
pData(dataBSmooth) <- pData

plotRegion(dataBSmooth,dmrs_sub[1,], extend = 5000, addRegions = dmrs_sub)

pdf(file = "output/bsmooth/dmrs_top_200_start_ordered.pdf", width = 10, height = 5)
plotManyRegions(dataBSmooth, dmrs_sub_ordered[1:200,], extend = 5000, addRegions = dmrs_sub_ordered)
dev.off()

save.image(file = "output/bsmooth/workspace.RData")

##############################################################
setwd("input_dir")

sim_chg_c20 <- read.delim("input_dir_for_input_files/C20p1_118_CHG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_t20 <- read.delim("input_dir_for_input_files/T20p1_117_CHG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_c20 <-  read.delim("input_dir_for_input_files/C20p1_118_CHH_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_t20 <-  read.delim("input_dir_for_input_files/T20p1_117_CHH_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_cpg_c20 <- read.delim("input_dir_for_input_files/C20p1_118_CpG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_cpg_c25 <- read.delim("input_dir_for_input_files/C25p2_126_CpG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))

dataBS_chg_c20 <- fitData(sim_chg_c20,"CHG")
dataBS_chg_t20 <- fitData(sim_chg_t20,"CHG")
dataBS_chh_c20 <- fitData(sim_chh_c20,"CHH")
dataBS_chh_t20 <- fitData(sim_chh_t20,"CHH")
dataBS_cpg_c20<-  fitData(sim_cpg_c20,"CpG")
dataBS_cpg_c25 <- fitData(sim_cpg_c25,"CpG")

peakRAM(DMR_calling(dataBS_chg_c20,"C20", "CHG")
peakRAM(DMR_calling(dataBS_chg_t20,"T20", "CHG")
peakRAM(DMR_calling(dataBS_chh_c20,"C20", "CHH")
peakRAM(DMR_calling(dataBS_chh_t20,"T20", "CHH")
peakRAM(DMR_calling(dataBS_cpg_c20,"C20", "CpG")
peakRAM(DMR_calling(dataBS_cpg_c25,"C25", "CpG")

