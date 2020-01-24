if (!require(methylSig)){
	if(!require(devtools)) install.packages(devtools)
	library(devtools)
	install_github('sartorlab/methylSig')
}

library(methylSig)
setwd("input_dir")
DMRcaller <- function(m_ecotype, m_context){
	newFolder <- paste(m_ecotype,"_",m_context,sep="")
	if (!dir.exists(paste("output/methylsig/",newFolder,sep=""))){
		dir.create(paste("output/methylsig/",newFolder,sep=""))
	}

	filelist <- list.files(paste("Data/methylkit/",newFolder,sep=""),
						pattern=paste(newFolder,"_methylKit_replicate","*", sep=""),
						full.names=T)

	sample_ids <- c(paste("test",seq(1:(length(filelist)/2)),sep=""),
					paste("control",seq(1:(length(filelist)/2)),sep=""))

	meth_cpg <- methylSigReadData(filelist, sample.ids = sample_ids, 
								assembly = "sim_cpg_g129", treatment = rep(c(1,0), 
								each = length(sample_ids)/2), minCount=7, context = "CpG",
								destranded = F)

	meth_cpg_tile <- methylSigTile(meth_cpg, win.size = 200)

	meth_cpg_dmr <- methylSigCalc(meth_cpg_tile, groups = c(1,0), min.per.group=3, 
								local.disp = T, winsize.disp= 200, 
							 	local.meth = T, winsize.meth = 200)

	meth_cpg_dmr_filtered <- meth_cpg_dmr[meth_cpg_dmr[,"qvalue"] < 0.05 & abs(meth_cpg_dmr[,"meth.diff"] > 10)]

	write.table(data.frame(meth_cpg_dmr_filtered@data.chr,
						meth_cpg_dmr_filtered@data.start,
						meth_cpg_dmr_filtered@data.end),
						paste("output/methylsig/",newFolder,"/methylsig_DMRs_",newFolder,"_200bp.bed",sep=""), 
						row.names=F,col.names=F,sep="\t",quote=F,dec = ",")
}
mergeDmrs <- function(m_ecotype,m_context){
	dmrs_to_merge <- read.table(paste("output/methylsig/",m_ecotype,"_",m_context,"/DMRs_200bp.bed",sep="")
								,colClasses=c("NULL","integer","integer"))
	colnames(dmrs_to_merge) <- c("start","end")
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
	mergedDMRs <- mergedDMRs[-1,]
	mergedDMRs <- cbind("chr1",mergedDMRs)


	write.table(mergedDMRs,paste("output/methylsig/",m_ecotype,"_",m_context,
								"/mergedDMRs_",m_ecotype,"_",m_context,sep=""),
				row.names=F,col.names=F,sep="\t",quote=F)
}

peakRAM(DMRcaller("Plus65", "CHH"))
peakRAM(DMRcaller("Minus65", "CHH"))
peakRAM(DMRcaller("Plus65", "CHG"))
peakRAM(DMRcaller("Minus65", "CHG"))
peakRAM(DMRcaller("Plus65", "CpG"))
peakRAM(DMRcaller("Minus65", "CpG"))

mergeDmrs("Plus65", "CpG")
mergeDmrs("Minus65", "CpG")
mergeDmrs("Plus65", "CHG")
mergeDmrs("Minus65", "CHG")
mergeDmrs("Plus65", "CHH")
mergeDmrs("Minus65", "CHH")

