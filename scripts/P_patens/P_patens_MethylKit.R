
if (!require(methylKit)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("methylKit")
}

library(methylKit)

fitData <- function(dfMethData,m_filename,m_context,dirName){
	methDf <- cbind(paste("chr1", dfMethData[,1],sep="_"),
					  rep("chr1",dim(dfMethData)[1]),
					  rep("+",dim(dfMethData)[1]),
					  dfMethData,stringsAsFactors=F)

	methDf <- methDf[,c(1,2,4,3,5:length(methDf))]

	newFolder <- file.path(paste(dirName,m_filename,"_",m_context,"/",sep=""))
	if (!dir.exists(newFolder)) {
		dir.create(newFolder)
	}
	for (i in 1:((length(methDf)-4)/2)) {
	write.table(cbind(methDf[c(1,2,3,4,i*2+3,i*2+4)], 100-methDf[i*2+4]), #adding also the absolute Thymine count
				paste(newFolder,m_filename,"_",m_context,"_","methylKit_replicate", i ,".txt", sep=""),
				col.names=c("id","chr","base","strand", "cov", "numC", "numT"),
				quote=F,
				row.names=F,
				sep="\t")
	}
}
DMR_caller <- function(m_ecotype, m_context){
	newFolder <- paste(m_ecotype,"_",m_context,"/",sep="")
	if (!dir.exists(paste("/homes/epi/nilay/1_Benchmarking/5_Output/Simulated/P_patens/methylkit",newFolder,sep=""))){
		dir.create(paste("/homes/epi/nilay/1_Benchmarking/5_Output/Simulated/P_patens/methylkit",newFolder,sep=""))
	}

	filelist <- list.files(paste("/data/MethylKit",m_ecotype,"_",m_context,sep=""),
						pattern=paste(m_ecotype,"_",m_context,"_methylKit_replicate*",sep=""),full.names=T)
	sample_list <- list()
	loc_list <- list()

	len_filelist <- length(filelist)

	for (i in 1:len_filelist){
		if (i <= (len_filelist/2)){
			loc_list[[i]] <- filelist[i]
			sample_list[[i]] <- paste("test",i, sep="")
		
		}else{
			loc_list[[i]] <- filelist[i]
			sample_list[[i]] <- paste("control",i - len_filelist/2, sep="")
		}
	}

	MethylRawList <- methRead( loc_list, sample.id=sample_list,
								pipeline=list(fraction=FALSE,chr.col=2,start.col=3,end.col=3,
								coverage.col=5,strand.col=4,freqC.col=6),
								treatment=rep(c(1,0),each=(len_filelist/2)),
								assembly="sim_plant_meth",context=m_context
								)

	MethylRawList_filtered <- filterByCoverage(MethylRawList, lo.count=10, hi.perc=99.9)
	MethylRawList_tiled <- tileMethylCounts(MethylRawList_filtered, win.size= 200, step.size= 200)
	MethylRawList_united <- unite(MethylRawList_tiled)
	MethylRawList_methDiffDSS <- calculateDiffMethDSS(MethylRawList_united, mc.cores=2)
	MethylRawList_methDiff_subsetDSS <- getMethylDiff(MethylRawList_methDiffDSS,difference=10,qvalue=0.05)

	sigDMR <- cbind(getData(MethylRawList_united[
							rownames(MethylRawList_united)%in%rownames(MethylRawList_methDiff_subsetDSS),]),
							getData(MethylRawList_methDiff_subsetDSS)
					)

	## Merge the DMRs
	dmrs_to_merge <- getData(MethylRawList_methDiff_subsetDSS)[2:3]
	mergedDMRs <- data.frame(start=0,end=0)
	insideDmr = 0
	for (i in 1:dim(dmrs_to_merge)[1]){
		if (i+1 > dim(dmrs_to_merge)[1]){
			if (insideDmr==1){
				mergedDMRs[dim(mergedDMRs)[1],2] <- dmrs_to_merge[i,2]
			}
			if (insideDmr==0){
				mergedDMRs <- rbind(mergedDMRs,dmrs_to_merge[i,1:2])
			}
		} else{
			if (dmrs_to_merge[i,2] == (dmrs_to_merge[i+1,1]-1)&& insideDmr==0) {
				mergedDMRs <- rbind(mergedDMRs,c(dmrs_to_merge[i,1],0))
				insideDmr = 1
			}
			if (dmrs_to_merge[i,2] != (dmrs_to_merge[i+1,1]-1)){
				
				if (insideDmr==0){
					mergedDMRs <- rbind(mergedDMRs,dmrs_to_merge[i,1:2])
				}
				if (insideDmr==1){
					mergedDMRs[dim(mergedDMRs)[1],2] <- dmrs_to_merge[i,2]
					insideDmr = 0
				}
			}
		}
	}

	# read out dmrs with additional information
	sigDMR <- sigDMR[-c((length(sigDMR)-6):(length(sigDMR)-3))]
	write.table(sigDMR,file=paste("/homes/epi/nilay/1_Benchmarking/5_Output/Simulated/P_patens/methylkit",newFolder,"DMR_",m_context,"_methylkitDSS.txt",sep=""),
				row.names=F,col.names=T,sep="\t",quote=F)

	options("scipen"=50) # avoid scientific name convention

	# read out the original dmrs (not merged)
	write.table(sigDMR[1:3],file = paste("/homes/epi/nilay/1_Benchmarking/5_Output/Simulated/P_patens/methylkit",newFolder,m_ecotype, "_DMRs_", m_context,".bed",sep=""),
				row.names=F,col.names=F,sep="\t",quote=F)

	# read out the merged dmrs
	mergedDMRs <- mergedDMRs[-1,]
	mergedDMRs <- cbind("chr1",mergedDMRs)
	write.table(mergedDMRs,file=paste("/homes/epi/nilay/1_Benchmarking/5_Output/Simulated/P_patens/methylkit",newFolder,"methylkit_DMRs_merged_",m_ecotype,"_",m_context,".bed",sep=""),
				row.names=F,col.names=F,sep="\t",quote=F)

	return(MethylRawList_methDiffDSS)
}


###########################################################
setwd("/input_dir")
library(methylKit)

sim_cpg_g129 <- read.delim("/data/G129_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_cpg_r143 <- read.delim("/data/R143_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_r143 <- read.delim("/data/R143_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_g129 <- read.delim("/data/G129_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_r143 <- read.delim("/data/R143_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_g129 <- read.delim("/data/G129_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))



fitData(sim_cpg_g129,"G129","CpG","/data/MethylKit")
fitData(sim_cpg_r143,"R143","CpG","/data/MethylKit")
fitData(sim_chg_r143,"R143","CHG","/data/MethylKit")
fitData(sim_chg_g129,"G129","CHG","/data/MethylKit")
fitData(sim_chh_r143,"R143","CHH","/data/MethylKit")
fitData(sim_chh_g129,"G129","CHH","/data/MethylKit")


g129_cpg_methdiffDSS<- peakRAM(DMR_caller("G129", "CpG"))
r143_cpg_methdiffDSS<- peakRAM(DMR_caller("R143", "CpG"))
r143_chg_methdiffDSS<- peakRAM(DMR_caller("R143", "CHG"))
g129_chg_methdiffDSS<- peakRAM(DMR_caller("G129", "CHG"))
r143_chh_methdiffDSS<- peakRAM(DMR_caller("R143", "CHH"))
g129_chh_methdiffDSS<- peakRAM(DMR_caller("G129", "CHH"))
