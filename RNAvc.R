script_name = "RNAvc.R"

args <- commandArgs(trailingOnly = FALSE)
args <- args[-1]


for (i in length(args):1) {
	if (any(grep(script_name, args[i]))) {
		root_dir <- args[i]
		root_dir <- sub(".*?/","/",root_dir,perl=TRUE)
		root_dir <- sub(script_name,"/",root_dir,perl=TRUE)
		args <- args[-i]
	}
	if (any(grep("CONFIG=", args[i]))) {
		config <- sub("CONFIG=","",args[i])
		args <- args[-i]
	}
}

config <- read.table(config)[,1]
threads <- config[pmatch("THREADS=",config)]
if (!is.null(threads) & !(is.na(threads))) {
	threads <- sub("THREADS=","",threads)
	threads=as.numeric(threads)
}

source(paste(root_dir,"RNAvc_lib.R",sep=""))

print (args)
if (length(args) > 0) {
	mode <- "train"
	features <- args[1]
} else {
	stop(paste("Usage:",args[1],args[2],"<feature file> [ <vcf file>  <model> ]"))
} 

if (length(args) > 1) {
	mode <- "predict"
	vcf <- args[2]
	model <- args[3]
}
save.image()

df.dirty <- read.table(features,sep="\t",header=TRUE)
#df <- cleanFeatures(df.dirty)
header <- sub(".features","",features)


#save.image()

if (mode == "train") {
	#model.ai <- paste(header,".atoi.model",sep="")
	#model.cu <- paste(header,".ctou.model",sep="")
	model.name <- paste(header,".model",sep="")

	df <- cleanFeatures(df.dirty)
	df2 <- df
	if (sum(is.na(df2$y)) > 0) { 
		df2[is.na(df2$y),]$y <- FALSE
		df2[df2$SiteClass == "Unknown",]$SiteClass <- "Edit"
	}

	#load("NA12878_big.atoi.model")

	gbm.error.model <- train.RNAvc(df2,cores=threads,class="error")
#save.image()

	# Calculate correctly classified errors
	output <- capture.output(test.RNAvc(gbm.error.model,df2,class="error"))
	cat(output,file=perf,sep="\n")
	# Calculate mis-classified edit sites
	if ("Edit" %in% names(table(df$SiteClass))) {
		output <- capture.output(test.RNAvc(gbm.error.model,df,class="AI"))
		cat(output,file=perf,sep="\n",append=TRUE)
	}

	df2.pred <- predict.RNAvc(gbm.error.model,df2,class="error")
	df3 <- df2.pred
	df3 <- df3[df3$RNAvcIsGood == TRUE,]
   	if ("RNAvcIsGood" %in% colnames(df3)) {
                df3 <- df3[-which(colnames(df3) == "RNAvcIsGood")]
        }
   	if ("RNAvcQualityScore" %in% colnames(df3)) {
                df3 <- df3[-which(colnames(df3) == "RNAvcQualityScore")]
        }


	gbm.ai.model <- train.RNAvc(df3,cores=threads)
	save(gbm.ai.model,gbm.error.model,file=model.name)
	#save(gbm.ai.model,file=model.ai)

#save.image()
	gbm.cu.model <- train.RNAvc(df3,cores=threads,class="CU")
	save(gbm.ai.model,gbm.cu.model,gbm.error.model,file=model.name)
save.image()

	if (require("caret") & require("e1071")) {
		perf <- paste(header,".perf",sep="")
		output <- capture.output(test.RNAvc(gbm.ai.model,df3,class="AI"))
		cat(output,file=perf,sep="\n")
		output <- capture.output(test.RNAvc(gbm.cu.model,df3,class="CU"))
		cat(output,file=perf,sep="\n",append=TRUE)
	}
} else if (mode == "predict") {
	df <- df.dirty
	vcf.annotated <- sub(".vcf",".annotated.vcf",vcf)
	load(model)
	predicted <- predict.models(df,error=gbm.error.model,ai=gbm.ai.model,cu=gbm.cu.model)
	#predicted <- predict.models(df,error=gbm.error.model)

	# Transitions are more likely to be true variants
	#predicted[(predicted$RNAvcQualityScore < 0.75) & (predicted$tstv == "TS"),]$RNAvcIsError <- FALSE
	#predicted <- predict.RNAvc(gbm.ai.model,predicted,class="AI")
	#predicted <- predict.RNAvc(gbm.cu.model,predicted,class="CU")
	annotateVCF(predicted,infile=vcf,outfile=vcf.annotated)
	if (length(table(df$SiteClass)) > 1) {
		perf <- paste(header,".perf",sep="")
		output <- capture.output(test.RNAvc(gbm.ai.model,predicted,class="AI"))
		cat(output,file=perf,sep="\n")
		predicted.cu <- predicted
		#predicted.cu[predicted.cu$SiteClass == "Unknown",]$y = TRUE
		predicted.cu[predicted.cu$SiteClass == "Unknown",]$SiteClass = "Edit"
		predicted.cu <- predicted.cu[!is.na(predicted.cu$y),]
		output <- capture.output(test.RNAvc(gbm.cu.model,predicted.cu,class="CU"))
		cat(output,file=perf,sep="\n",append=TRUE)
		#output <- capture.output(test.RNAvc(gbm,predicted))
		#cat(output,perf,sep="\n")
	}
}


