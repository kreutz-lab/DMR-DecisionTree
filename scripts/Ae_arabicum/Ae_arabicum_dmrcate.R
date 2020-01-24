#14.11.17 

if (!require(DMRcate)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("DMRcate")
}


fitData <- function(dfMethData,m_context){
	methDf <- cbind(c(rep("chr1",dim(dfMethData)[1])),
					c(rep("+",dim(dfMethData)[1])),
					c(rep(m_context,dim(dfMethData)[1])),
					dfMethData,stringsAsFactors=F)

	colnames(methDf)[1:4] <- c("chr","strand","context","pos")

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

