if (!require(DMRcate)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("DMRcate")
}


fitData <- function(dfMethData,m_context){
	## In order to read the data, the raw data have to be adjusted
	# chromosome name, strand and context is missing
	methDf <- cbind(c(rep("chr1",dim(dfMethData)[1])),
					c(rep("+",dim(dfMethData)[1])),
					c(rep(m_context,dim(dfMethData)[1])),
					dfMethData,stringsAsFactors=F)

	colnames(methDf)[1:4] <- c("chr","strand","context","pos")

	#change the relation to absolute values
	for (i in 1:(length(methDf)/2-2)){
		methDf[i*2+4] <- methDf[i*2+3]*((methDf[i*2+4]/100))
		if(i<2){
			colnames(methDf)[(i*2+3):(i*2+4)] <- c(paste("Treatment",i,".cov",sep=""),paste("Treatment",i,".C",sep=""))
		}else{
			colnames(methDf)[(i*2+3):(i*2+4)] <- c(paste("Control",i-3,".cov",sep=""),paste("Control",i-3,".C",sep=""))

		}
	}

	# Creates a Grange class
	BS <- makeGRangesFromDataFrame(methDf, seqnames.field="chr", 
						start.field="pos", 
						end.field="pos", 
						strand.field="strand", 
						keep.extra.columns=T)
	return(BS)
}

DMR_calling <- function(m_dataMeth,m_ecotype, m_context){
	newFolder <- paste(m_ecotype,"_",m_context,sep="")
	if (!dir.exists(paste("output/dmrcate/",newFolder,sep=""))){
		dir.create(paste("output/dmrcate/",newFolder,sep=""))
	}

	coverage <- as.data.frame(m_dataMeth)[,c(1:2, grep(".cov$", colnames(as.data.frame(m_dataMeth))))]
	meth <- as.data.frame(m_dataMeth)[,c(1:2, grep(".C$", colnames(as.data.frame(m_dataMeth))))]

	treat1 <- data.frame(chr= coverage$seqnames, pos=coverage$start, N=coverage$Treatment1.cov, X=meth$Treatment1.C)
	treat2 <- data.frame(chr= coverage$seqnames, pos=coverage$start, N=coverage$Treatment2.cov, X=meth$Treatment2.C)
	treat3 <- data.frame(chr= coverage$seqnames, pos=coverage$start, N=coverage$Treatment3.cov, X=meth$Treatment3.C)

	ctrl1 <- data.frame(chr= coverage$seqnames, pos=coverage$start, N=coverage$Control1.cov, X=meth$Control1.C)
	ctrl2 <- data.frame(chr= coverage$seqnames, pos=coverage$start, N=coverage$Control2.cov, X=meth$Control2.C)
	ctrl3 <- data.frame(chr= coverage$seqnames, pos=coverage$start, N=coverage$Control3.cov, X=meth$Control3.C)

	samples <- list(treat1,treat2,treat3,ctrl1,ctrl2,ctrl3)
	sampNames <- sub("\\..*","",colnames(meth))[-c(1:2)]

	obj_bsseq <- makeBSseqData(samples, sampNames)

	DSSres <- DMLtest(obj_bsseq, group1=sampNames[1:3], group2=sampNames[4:6], smoothing= F)

	wgbsannot <- cpg.annotate(datatype="sequencing", DSSres,analysis.type="differential")
	# DMR calling:
	wgbs.DMRs <- dmrcate(wgbsannot, lambda = 1000, C = 50, pcutoff = 0.05, mc.cores = 1)

	## extract DMR positions and then save the data
	df = transform(wgbs.DMRs$results,coord=colsplit(coord, split = "\\-", names= c("a","b")))

	startDmr <- c()
	for (i in 1:dim(df$coord)[1]){
		startDmr[i] <- as.integer(substr(as.character(df$coord$a[i]),6,nchar(as.character(df$coord$a[i]))))
	}

	endDmr <- as.vector(df$coord$b)
	data_dmr_to_bed <- data.frame(chr=rep("chr1",length(startDmr)),startDmr,endDmr,stringsAsFactors=F)
	data_dmr_to_bed <- data_dmr_to_bed[order(as.integer(data_dmr_to_bed$startDmr)),]

	write.table(data_dmr_to_bed,paste("output/dmrcate/",newFolder,"/DMRcate_dmrs_",m_ecotype,"_",m_context,".bed",sep=""),
				row.names=F,col.names=F,sep="\t",quote=F)
}

# plot DMR with DMRcate plotting methods
wgbs.ranges <- extractRanges(wgbs.DMRs, genome= "simulated")
groups <- c(Treatment = "darkorange", Control = "blue")
cols <- groups[sub("[0-9]","",sampNames)]
DMR.plot(ranges = wgbs.ranges, dmr = 1, CpGs = m_dataMeth, phen.col = cols, genome = "simulated")

###############################################
setwd("input_dir")

library(DMRcate)
library(reshape)


sim_cpg_g129 <- read.delim("Data/G129_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),2)))
sim_cpg_r143 <- read.delim("Data/R143_CpG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_r143 <- read.delim("Data/R143_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chg_g129 <- read.delim("Data/G129_CHG_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_r143 <- read.delim("Data/R143_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))
sim_chh_g129 <- read.delim("Data/G129_CHH_SimulationData.txt",colClasses=c("integer",rep(c("integer","numeric"),6)))

dataMeth_cpg_g129 <- fitData(sim_cpg_g129,"cpg")
dataMeth_cpg_r143 <- fitData(sim_cpg_r143,"cpg")
dataMeth_chg_r143 <- fitData(sim_chg_r143,"chg")
dataMeth_chg_g129 <- fitData(sim_chg_g129,"chg")
dataMeth_chh_r143 <- fitData(sim_chh_r143,"chh")
dataMeth_chh_g129 <- fitData(sim_chh_g129,"chh")

peakRAM(DMR_calling(dataMeth_cpg_g129,"G129", "CpG"))
peakRAM(DMR_calling(dataMeth_cpg_r143,"R143", "CpG"))
peakRAM(DMR_calling(dataMeth_chg_r143,"R143", "CHG"))
peakRAM(DMR_calling(dataMeth_chg_g129,"G129", "CHG"))
peakRAM(DMR_calling(dataMeth_chh_r143,"R143", "CHH"))
peakRAM(DMR_calling(dataMeth_chh_g129,"G129", "CHH"))

