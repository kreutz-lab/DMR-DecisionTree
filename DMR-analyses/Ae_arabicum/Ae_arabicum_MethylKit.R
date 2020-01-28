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
	if (!dir.exists(paste("input_dir",newFolder,sep=""))){
		dir.create(paste("input_dir",newFolder,sep=""))
	}

	filelist <- list.files(paste("input_dir/",m_ecotype,"_",m_context,sep=""),
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
	MethylRawList_united <- unite(MethylRawList_tiled,min.per.group = 1L)
	MethylRawList_methDiffDSS <- calculateDiffMeth(MethylRawList_united, mc.cores=2)
	MethylRawList_methDiff_subsetDSS <- getMethylDiff(MethylRawList_methDiffDSS,difference=10,qvalue=0.05)



	sigDMR <- cbind(getData(MethylRawList_united[
							rownames(MethylRawList_united)%in%rownames(MethylRawList_methDiff_subsetDSS),]),
							getData(MethylRawList_methDiff_subsetDSS)
					)

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

	sigDMR <- sigDMR[-c((length(sigDMR)-6):(length(sigDMR)-3))]
	write.table(sigDMR,file=paste("input_dir",newFolder,"DMR_",m_context,"_methylkitDSS.txt",sep=""),
				row.names=F,col.names=T,sep="\t",quote=F)

	options("scipen"=50) # avoid scientific name convention

	write.table(sigDMR[1:3],file = paste("input_dir",newFolder,m_ecotype, "_DMRs_", m_context,".bed",sep=""),
				row.names=F,col.names=F,sep="\t",quote=F)

	# read out the merged dmrs
	mergedDMRs <- mergedDMRs[-1,]
	mergedDMRs <- cbind("chr1",mergedDMRs)
	write.table(mergedDMRs,file=paste("input_dir",newFolder,"methylkit_DMRs_merged_",m_ecotype,"_",m_context,".bed",sep=""),
				row.names=F,col.names=F,sep="\t",quote=F)

	return(MethylRawList_methDiffDSS)
}


###########################################################
setwd("input_dir_for_input_files/")
library(methylKit)
sim_chg_c20 <- read.delim("input_dir_for_input_files/C20p1_118_CHG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_t20 <- read.delim("input_dir_for_input_files/T20p1_117_CHG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_c20 <-  read.delim("input_dir_for_input_files/C20p1_118_CHH_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_t20 <-  read.delim("input_dir_for_input_files/T20p1_117_CHH_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_cpg_c20 <- read.delim("input_dir_for_input_files/C20p1_118_CpG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_cpg_c25 <- read.delim("input_dir_for_input_files/C25p2_126_CpG_Sc65_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))


fitData(sim_chg_c20p,"C20","CHG","input_dir")
fitData(sim_chg_t20p,"T20","CHG","input_dir")
fitData(sim_chh_t20,"T20","CHH","input_dir")
fitData(sim_chh_c20,"C20","CHH","input_dir")
fitData(sim_cpg_c20p,"C20","CpG","input_dir")
fitData(sim_cpg_c25p,"C25","CpG","input_dir")

c20_chg_methdiffDSS=peakRAM(DMR_caller"C20","CHG"))
t20_chg_methdiffDSS=peakRAM(DMR_caller"T20","CHG"))
t20_chh_methdiffDSS=peakRAM(DMR_caller"T20","CHH"))
c20p_chh_methdiffDSS=peakRAM(DMR_caller"C20","CHH"))
c20p_cpg_methdiffDSS=peakRAM(DMR_caller"C20","CpG"))
c25p_cpg_methdiffDSS=peakRAM(DMR_caller"C25","CpG"))


