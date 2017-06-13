LOCALTEST <- FALSE

if (LOCALTEST) {
	library(qqman)
	setwd("/home/mbrown/DNAnexus/")
	source("mhplot/resources/mhplot_functions.R")	

	# Testing variables
	data.files <- Sys.glob("data/mb_test_data_chr*.Rda")
	# snpinfo.files <- Sys.glob("/home/klawson/WGS_frz3/snpinfo/nomonom/EA*snpinfo*chr22*RData")
	p <- "p"
	gene <- "gene"
	chr <- "chr"
	pos <- "pos"
	snp <- "Name"
	filter <- "p<=0.01"
	mh.color <- c("gray50","black")
	chrlabs <- NULL
	suggestiveline <- -log10(1e-05)
	genomewideline <- -log10(5e-08)
	highlight <- NULL
	logp <- TRUE
	main <- "Main Title Here"
	output.type <- "png"

} else {
	## Parse the arguments
	args <- commandArgs(trailingOnly=TRUE)
	args <- paste(args,collapse=" ")
	print(args)
	args <- unlist(strsplit(args,";"))
	print(args)

	## Mandatory parameters
	data.files <- unlist(strsplit(args[14]," "))

	## Optional parameters
	# snpinfo.files <- 
	p <- args[1]
        gene <- args[2]
        chr <- args[3]
        pos <- args[4]
        snp <- args[5]
        chrlabs <- args[6]
        # offset <- as.numeric(args[7])
	filter <- args[7]
        suggestiveline <- eval(parse(text=args[8]))
        genomewideline <- eval(parse(text=args[9]))
        gene.list <- unlist(strsplit(args[10]," "))
        logp <- args[11]
        main <- args[12]
        mh.color <- unlist(strsplit(args[13]," "))
	output.type <- args[15]

	## Install libraries
	install.packages("/qqman_0.1.2.tar.gz", repos=NULL, type="source")
	library(qqman)
}

# Functions
get.delim <- function(f) {
        d <- readLines(f,n=1)
        delim <- ifelse(grepl(" ", d)," ",ifelse(grepl(",",d),",","\t"))
        return(delim)
}

## Load data
l <- vector(mode="list",length=length(data.files))
# Trim any whitespace in the file names
data.files <- sapply(data.files, function(x) gsub("^\\s+|\\s+$","",x))
for (idx in 1:length(data.files)) {
	cat("Loading ",data.files[idx],"\n")
	if (grepl("Rda$", data.files[idx])) {
		load(data.files[idx])
		rm(skat)
	} else {
		delim <- get.delim(data.files[idx])
		sing <- read.table(data.files[idx],header=T,as.is=T,sep=delim)
	}
	
	if (idx==1) {
		print(head(sing))
	}

	## Use filter string to subset the data
	sing <- eval(parse(text = paste0("subset(sing,",filter,")")))

	## Check other variables
	chrlabs <- NULL

	## Add chr and pos if not in the dataset and can be extracted from snp name
	if (!is.element(chr,names(sing))) {
		if (grepl(":",sing[1,snp])) {
			message("Extracting chr from SNP name...")
			sing[,chr] <- as.numeric(gsub(":[0-9]+$","",sing[,snp]))
		} else {
			message("Extracting pos from SNP name...")
			stop("Dataset does not contain valid chr column and cannot be extracted from SNP name.")
		}
	}
	if (!is.element(pos,names(sing))) {
		if (grepl(":",sing[1,snp])) {
			sing[,pos] <- as.numeric(gsub("[0-9]+:","",sing[,snp]))
		} else {
			stop("Dataset does not contain valid position column and cannot be extracted from SNP name.")
		}
	}

	# Subset to save memory
	wh.col <- which(c(p,gene,chr,pos,snp) %in% names(sing))
	l[[idx]] <- sing[!is.na(sing[,p]),c(p,gene,chr,pos,snp)[wh.col]]
	rm(sing)
}
gc()

# Combine pieces into one dataset. Should be fine if data is all in one file
message("Merging chromosomes...")
d <- do.call("rbind",l)
rm(l)
gc()

## Variables
# gene.list <- c("chr9_win7788","chr12_win32617","chr16_win6755")
if (length(gene.list)==1) {
	if (gene.list == "NULL") {
		highlight <- NULL
	} else {
		if ( is.element(gene,names(d))) {
                	highlight <- d[d[,gene] %in% gene.list,"Name"]
        	} else {
                	highlight <- NULL
        	}
	}
} else {
	if ( is.element(gene,names(d))) {
		highlight <- d[d[,gene] %in% gene.list,"Name"]
	} else {
		highlight <- NULL
	}
}
	
## Plot the data
message("Plotting manhattan")
eval(parse(text=paste0(output.type,"(\"mh.",output.type,"\",height=480,width=720,units=\"px\")")))
# png("test_mh.png",height=10,width=16,units="cm",res=300)
# This was the data.table code
# manhattan(d[d$p<=offset,], chr=chr, bp=pos, p=p, snp=snp, col=mh.color, chrlabs=chrlabs,
#        suggestiveline=suggestiveline, genomewideline=genomewideline, highlight=highlight, logp=logp, main=main)
	
# print(c(chr=chr, bp=pos, p=p, snp=snp, col=mh.color, chrlabs=chrlabs,
#       suggestiveline=suggestiveline, genomewideline=genomewideline, highlight=highlight, logp=logp, main=main))
# dd <- na.omit(d[d$p<=offset,])
# manhattan(na.omit(d[d$p<=offset,]), chr=chr, bp=pos, p=p, snp=snp, col=mh.color, chrlabs=chrlabs,
manhattan(na.omit(d), chr=chr, bp=pos, p=p, snp=snp, col=mh.color, chrlabs=chrlabs,
	       suggestiveline=suggestiveline, genomewideline=genomewideline, highlight=highlight, logp=logp, main=main)
dev.off()
	



