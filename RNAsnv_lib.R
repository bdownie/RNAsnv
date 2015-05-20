
train.RNAsnv <- function (train,minimumCumulativeImportance=100,dist="bernoulli",save.model=NULL,cores=5,save.glm=NULL,class="AI",max.trees=50000) {

	#train <- cleanFeatures(train)
	#train <- train[!is.na(train$y),]
	#train <- train[train$MutationType == "AG" | train$MutationType == "TC",]
	if (class == "error") { 
		train <- train[get_tstv(train$MutationType) == "TV",]
	} else {
		train <- only_edit(train,class=class)
	}

	if (sum(train$SiteClass == "Unknown") > 0) { 
		train[train$SiteClass == "Unknown",]$y <- FALSE
	}
		
	if ("site.id" %in% colnames(train)) {
		train <- train[-which(colnames(train) == "site.id")]
	}
	if ("Sample" %in% colnames(train)) {
		train <- train[-which(colnames(train) == "Sample")]
	}
	if ("SiteClass" %in% colnames(train)) {
		train <- train[-which(colnames(train) == "SiteClass")]
	}
	if ("IndbSNP" %in% colnames(train)) {
		train <- train[-which(colnames(train) == "IndbSNP")]
	}
	if ("MutationType" %in% colnames(train)) {
		train <- train[-which(colnames(train) == "MutationType")]
	}


	# Train variants only on high variant quality


	ntree <- 400
	best.iter <- ntree

	train.bak <- train
	sample_count = 10000
	#set.seed(1)
	# Optimize for subset of features
	indices <- logical(length=ncol(train))

	require("gbm")

	while ((sum(indices) < length(indices)) || ((best.iter/ntree) > 0.8)) {
		if (ntree > (max.trees/5)) { 
			ntree = max.trees
		} else { 
			ntree = ntree * 5
		}

		#if (ntree > max.trees) { 
		#	ntree = max.trees
		#}

		if (nrow(train) > sample_count) { 
			train.sample <- train[sample(1:nrow(train),size=sample_count),]
		} else { train.sample = train }

		gbm.model <- gbm(y~.,data=train.sample,distribution=dist,n.trees=ntree,
			shrinkage=0.01, interaction.depth=5, bag.fraction = 0.5, train.fraction = 0.5,
			n.minobsinnode = 10, cv.folds = 5,keep.data=TRUE,verbose=FALSE,n.cores=cores)

		# Discard incidental features from each iteration (to speed up iteration computation)
		best.iter <- gbm.perf(gbm.model,method="cv",plot=FALSE)
		cumulative <- summary(gbm.model,n.trees=best.iter,plot=FALSE)
		cumulative <- cumulative[cumulative$rel.inf > 1,]
		indices <- logical(length=ncol(train))
		names(indices) <- colnames(train)

		for (i in 1:nrow(cumulative)) {
			indices[[rownames(cumulative)[i]]] <- TRUE
		}
		indices["y"] <- TRUE
		if (sum(indices) < length(indices)) { 
			train <- subset(train,select=indices)
		} else {
			if (ntree >= max.trees) { 
				best.iter = 0
			}
		#	 else if (ntree > (max.trees/5)) { 
		#		ntree = max.trees
		#	} else { 
		#		ntree = ntree * 5
		#	}
		}
		#train.bak <- subset(train.bak,select=indices)
	}

	ntree <- 400
	best.iter <- ntree
	while ((sum(indices) < length(indices)) || ((best.iter/ntree) > 0.8)) {
		if (ntree > (max.trees/5)) { 
			ntree = max.trees
		} else { 
			ntree = ntree * 5
		}

	#while ((best.iter/ntree) > 0.8) {
		#if (ntree > max.trees) { 
		#	ntree = max.trees
		#}

		train.sample <- train


		gbm.model <- gbm(y~.,data=train.sample,distribution=dist,n.trees=ntree,
			shrinkage=0.005, interaction.depth=5, bag.fraction = 0.5, train.fraction = 0.5,
			n.minobsinnode = 10, cv.folds = 5,keep.data=TRUE,verbose=FALSE,n.cores=cores)

		# Discard incidental features from each iteration (to speed up iteration computation)
		best.iter <- gbm.perf(gbm.model,method="cv",plot=FALSE)
		cumulative <- summary(gbm.model,n.trees=best.iter,plot=FALSE)
		cumulative <- cumulative[cumulative$rel.inf > 1,]
		indices <- logical(length=ncol(train))
		names(indices) <- colnames(train)

		for (i in 1:nrow(cumulative)) {
			indices[[rownames(cumulative)[i]]] <- TRUE
		}
		indices["y"] <- TRUE
		if (sum(indices) < length(indices)) { 
			train <- subset(train,select=indices)
		} else {
			if (ntree >= max.trees) { 
				best.iter = 0
			}
		}
	}
	#train <- train.bak

	#ntree = best.iter

	#while ((best.iter/ntree) > 0.8) {
	#	ntree = best.iter * 2 
	#	if (ntree > max.trees) { ntree = max.trees }
#
#		gbm.model <- gbm(y~.,data=train,distribution=dist,n.trees=ntree,
#			shrinkage=0.005, interaction.depth=5, bag.fraction = 0.5, train.fraction = 0.5,
#			n.minobsinnode = 10, cv.folds = 5,keep.data=TRUE,verbose=FALSE,n.cores=cores)
#		best.iter <- gbm.perf(gbm.model,method="cv",plot=FALSE)
#		if (ntree == max.trees) { best.iter = 0 }
#	}
	#if (!is.null(save.glm)) {
	#	glm.model <- glm(y~.,data=train,family=binomial(logit))
	#	save(glm.model,file=save.glm)
	#}
			

	if (!is.null(save.model)) { 
		save(gbm.model,file=save.model)
	}
	gbm.model

	#return_val
}

predict.RNAsnv <- function (gbm.model,df,ntree=NULL,class="AI",decision.boundry=0.5) {
	library(gbm)
	#rm(df.tmp)
	if (!("response" %in% colnames(df))) {
		df$response <- 1
	}

	if ("MutationType" %in% colnames(df)) {
		if (class == "AI") {
			index <-  ((df$MutationType == "AG") | (df$MutationType == "TC"))
		} else if (class == "CU") { 
			index <-  ((df$MutationType == "GA") | (df$MutationType == "CT"))
		} else { 
			index <- rep(TRUE,nrow(df))
		}		
		df[index,]$response <- predict(gbm.model,df,type="response")[index]
		#df[get_tstv(df$MutationType) == "TV",]$response = 0
		#df[!(df$MutationType == "AG" | df$MutationType == "TC"),]$response = 0

	}
	df$Predicted <- logical(length=nrow(df))
	#df[df$response >= 0.5,]$Predicted = TRUE
	if (class == "error") { 
		colnames(df)[which(colnames(df) == "Predicted")] <- "RNAsnvIsGood"
        	colnames(df)[which(colnames(df) == "response")] <- "RNAsnvQualityScore"
		df[get_tstv(df$MutationType) == "TS" & df$RNAsnvQualityScore >= 0.25,]$RNAsnvIsGood <- TRUE
		df[get_tstv(df$MutationType) == "TV" & df$RNAsnvQualityScore >= decision.boundry,]$RNAsnvIsGood <- TRUE
	} else {
		if (!("RNAsnvScore" %in% colnames(df))) {
			df$RNAsnvScore = 1
		}
		if (!("RNAsnvClass" %in% colnames(df))) {
			df$RNAsnvClass <- factor("SNV",levels=c("Edit","SNV"))
		}
		df$RNAsnvScore = pmin(df$RNAsnvScore,df$response)
		df[df$RNAsnvScore >= decision.boundry,]$RNAsnvClass <- "SNV"
		df[df$RNAsnvScore < decision.boundry,]$RNAsnvClass <- "Edit"
	}



	df
	#list(df,df.retain)
}



test.RNAsnv <- function (gbm.model,test, response=FALSE,class="AI") { 
	library(gbm)
	library(caret)

	test <- predict.RNAsnv(gbm.model,test,class=class)

	print(paste("Confusion matrix for edit class:",class))

	if (class == "error") {
		test2 <- test[get_tstv(test$MutationType) == "TV",]
		if (table(test2$SiteClass)["Unknown"] > 0) { 
			test2[test2$SiteClass == "Unknown",]$y <- FALSE
		}
		test2$SiteClass <- droplevels(test2$SiteClass)
		print(confusionMatrix(test2$RNAsnvIsGood,test2$y,positive="TRUE"))
	}
	else { 
		test2 <- only_edit(test,class=class,keep.MutationType=TRUE)
		test2$SiteClass <- droplevels(test2$SiteClass)
		print(confusionMatrix(test2$RNAsnvClass,test2$SiteClass,positive="SNV"))
	}

	if (response==TRUE) { 
		test
	}
}

cleanFeatures <- function (df,sample=NULL) {

	rownames(df) <- 1:nrow(df)


	# Normalize distances to top99% of length
	median_99quant <- quantile(df$MedianDistFromEnd,probs=0.99,na.rm=TRUE)
	max_99quant <- quantile(df$MaxDistFromEnd,probs=0.99,na.rm=TRUE)
	df$MedianDistFromEnd[df$MedianDistFromEnd > median_99quant] <- median_99quant
	df$MaxDistFromEnd[df$MaxDistFromEnd > max_99quant] <- max_99quant


	#df$RNAsnv_score <- (df$MedianDistFromEnd + df$MaxDistFromEnd)

	# Remove > 1 alternate locus
	df$MutationType[grep(",",df$MutationType)] <- NA
	df$MutationType <- factor(df$MutationType)

	# Need to remove microsatellites
	df <- df[df$RepeatClass != "Simple_repeat",]
	df <- df[df$RepeatClass != "Low_complexity",]
	df$RepeatClass <- droplevels(df$RepeatClass)

	if (!"Unknown" %in% levels(df$Biotype)) { 
		df$Biotype <- factor(df$Biotype,levels=c(levels(df$Biotype),"Unknown"))
	}
	df$Biotype[is.na(df$Biotype)] <- "Unknown"
	df$SpliceDist[is.na(df$SpliceDist)] <- -1
	df$RepeatDist[is.na(df$RepeatDist)] <- -1

	# Remove columns with only NA
	df <- df[,colSums(is.na(df)) != nrow(df)]
	#df <- df[!is.na(df$SiteClass),]
	if ("SiteClass" %in% colnames(df)) {
		df <- df[!is.na(df$SiteClass),]
	}
	df <- df[,colSums(is.na(df)) < nrow(df)/10]
	df <- df[rowSums(is.na(df)) == 0,]

	df$tstv <- get_tstv(df$MutationType)
	df$cv_tstv <- get_tstv(df$ClosestVarType)

	#df$is.edit <- getSiteClass(df,
	#factor("Unknown",levels=c("Unknown","SNP","Edit"))
	#df$is.edit[df$Verified == "Verified"] <- "Edit"
	#df$is.edit[df$DNASite == "Verified"] <- "SNP"

	df$site.id <- paste(df$Chr,df$Site,sep=".")

	# Clean all the columns which should go away
	if ("Chr" %in% colnames(df)) {
		df <- df[-which(colnames(df) == "Chr")]
	}
	if ("Site" %in% colnames(df)) {
		df <- df[-which(colnames(df) == "Site")]
	}
	if ("GeneID" %in% colnames(df)) {
		df <- df[-which(colnames(df) == "GeneID")]
	}

	cn <- colnames(df)
#	if ("site.id" %in% cn) {
#		cn <- cn[-which(cn == "site.id")]
#	}
	# Remove all columns with variance of 0.
	for (i in (length(cn) - 1):1) {
		if (is.numeric(df[,i]) & (var(df[,i],na.rm=TRUE) ==  0)) {
			df <- df[-i]
		} 
		else if (is.logical(df[,i])) {
			df[,i] <- as.factor(df[,i])
		}
		if (is.factor(df[,i]) & (length(table(df[,i],useNA="no")) == 1)) {
			df <- df[-i]
		}
	}
		# limit factors to no larger than 53 (for randomforest)
		# This may cause problems between test and train data set, so should handle filtering at feature extraction.
			#else if (is.factor(df[i])) {
			#x <- summary(df[i])
			#if (length(x) > max_factor_levels) {
			#	x <- x[rev(order(x))]
			#	delete_levels <- names(x[(max_factor_levels+1):length(x)])
			#	if ("Unknown" %in% names(x)) {
			#		del_val <- "Unknown"
			#	} else {
			#		del_val = NA
			#	}
			#	df[df[,i] %in% delete_levels,i] <- del_val
			#	df[,i] <- droplevels(df[,i])
			#}
		#}



#	df.retain <- df[df$SiteClass == "Unknown",]
	df$y <- logical(length=nrow(df))
	if ("SiteClass" %in% colnames(df)) {
		if (length(df[df$SiteClass == "Edit",]$SiteClass) > 0) {
			df[df$SiteClass == "SNV",]$y <- TRUE
		}
		if (sum(is.na(df$SiteClass)) > 0) { 
			df[df$SiteClass == "Unknown",]$y <- NA
		}
	}

	# Remove rows with at least one NA
	df 
	#list(df,df.retain,names=c("Cleaned","Retained"))
}

get_tstv <- function(v) { 
	output <- factor(rep("TV",length(v)),levels=c("TS","TV"))
	#output <- "TV"
	output[v == "AG"] <- "TS"
	output[v == "GA"] <- "TS"
	output[v == "CT"] <- "TS"
	output[v == "TC"] <- "TS"

	output
}


only_edit <- function(df,class="AI",keep.MutationType=FALSE) { 
	if (!("MutationType" %in% colnames(df))) {
		stop("MutationType column not found!")
	}
	if (class == "AI") { 
		df <- df[df$MutationType == "AG" | df$MutationType == "TC",]
	} else if (class == "CU") { 
		df <- df[df$MutationType == "GA" | df$MutationType == "CT",]
	#	if ("SiteClass" %in% colnames(df)) {
	#		df[df$SiteClass != "SNV",]$SiteClass <- "Edit"
	#	}
	} else { 
		stop("Unknown edit class")
	}

	if ("SiteClass" %in% colnames(df)) {
		if (table(df$SiteClass)["Edit"] == 0) {
			df[df$SiteClass != "SNV",]$SiteClass <- "Edit"
		}
		df <- df[df$SiteClass != "Unknown",]
	}

	if ("tstv" %in% colnames(df)) {
		df <- df[-which(colnames(df) == "tstv")]
	}
	if ("MutationType" %in% colnames(df) & !keep.MutationType) {
		df <- df[-which(colnames(df) == "MutationType")]
	}
	if ("y" %in% colnames(df)) {
		df[df$SiteClass == "SNV",]$y  <- TRUE
		df <- df[!is.na(df$y),]
	}

	df
}

annotateVCF <- function(prediction,infile=NULL,outfile=NULL) {
	header <- system(paste("grep ^#",infile),intern=TRUE)
	lastline <- tail(header,n=1)

	header=head(header,n=(length(header) - 1))
	RNAsnv_text = c(
		"##RNAsnvVersion=1.0",
#		"##INFO<ID=MedianEditDist,Number=1,Type=Integer,Description=\"Median edit distance as reported by STAR for all reads at site.\">",
#		"##INFO<ID=MedianReadEndDist,Number=1,Type=Integer,Description=\"Median distance of variant from the nearest read end.\">",
#		"##INFO<ID=SNPClust,Number=1,Type=Integer,Description=\"Number of SNPs in 10 bp window.\">",
		"##INFO<ID=RNAsnvIsGood,Number=1,Type=Character,Description=\"Classification of variant as error.\">",
		"##INFO<ID=RNAsnvQualityScore,Number=1,Type=Character,Description=\"Confidence of classification of variant as an error (0 = SNV, 1=Edit site).\">",
		"##INFO<ID=RNAsnvClass,Number=1,Type=Character,Description=\"Classification of variant as edit site, genomic SNV, or Not predicted.\">",
		"##INFO<ID=RNAsnvScore,Number=1,Type=Character,Description=\"Confidence of classification (0 = SNV, 1=Edit site).\">"
	)

	header <- c(header,RNAsnv_text,lastline)

	df <- read.table(infile)
	tmp <- paste(df$V1,df$V2,sep=".")
	df <- df[!duplicated(tmp),]
	rownames(df) <- tmp[!duplicated(tmp)]

	RNAsnv <- data.frame( RNAsnvQualityScore=rep(NA,nrow(df)), RNAsnvIsGood=factor(rep("NotPredicted",nrow(df)),levels=c("NotPredicted","TRUE","FALSE")), RNAsnvScore=rep(NA,nrow(df)), RNAsnvClass=factor(rep("NotPredicted",nrow(df)),levels=c("NotPredicted","Edit","SNV")))
	rownames(RNAsnv) <- rownames(df)


	if (!("site.id" %in% colnames(prediction))) { 
		if (("Chr" %in% colnames(prediction)) & ("Site" %in% colnames(prediction))) { 
			prediction$site.id <- paste(prediction$Chr,prediction$Site,sep=".")
		} else {
			stop("Couldn't find or generate site.id column for annotating VCF")
		}
	}
	prediction <- prediction[!duplicated(prediction$site.id),]
	rownames(prediction) <- prediction$site.id


	index <- pmatch(rownames(prediction),rownames(RNAsnv))
	RNAsnv[index,1] <- signif(prediction$RNAsnvQualityScore,digits=3)
	RNAsnv[index,2] <- prediction$RNAsnvIsGood	
	RNAsnv[index,3] <- signif(prediction$RNAsnvScore,digits=3)
	RNAsnv[index,4] <- prediction$RNAsnvClass

	#append.score[1:length(append.score)] <- signif(prediction[index,]$response,digits=3)

	#append.class <- factor(rep("NotPredicted",nrow(df)),levels=c("NotPredicted","EditSite","SNV"))
	#names(append.class) <- rownames(df)

	#trues <- prediction[prediction$Predicted == TRUE,]$site.id
	#falses <- prediction[prediction$Predicted != TRUE,]$site.id
	#append.class[names(append.class) %in% trues] <- "EditSite"
	#append.class[names(append.class) %in% falses] <- "SNV"
#
#	append.score <- as.numeric(rep(NA,length=nrow(df)))
#	names(append.score) <- rownames(df)


#	rownames(prediction) <- prediction$site.id
#	index <- pmatch(names(append.score),rownames(prediction))
	#append.score[1:length(append.score)] <- signif(prediction[index,]$response,digits=3)

	df[,8] <- paste(df[,8],";RNAsnvIsGood=",RNAsnv$RNAsnvIsGood,sep="")
	df[,8] <- paste(df[,8],";RNAsnvQualityScore=",RNAsnv$RNAsnvQualityScore,sep="")
	df[,8] <- paste(df[,8],";RNAsnvClass=",RNAsnv$RNAsnvClass,sep="")
	df[,8] <- paste(df[,8],";RNAsnvScore=",RNAsnv$RNAsnvScore,sep="")
	
	write.table(header,file=outfile,append=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	write.table(df,file=outfile,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}


predict.models <- function(df,error=NULL,ai=NULL,cu=NULL) {
	if (!is.null(error)) {
		df <- predict.RNAsnv(error,df,class="error")
	}
	if (!is.null(ai)) {
		df<- predict.RNAsnv(ai,df,class="AI")
	}
	if (!is.null(cu)) {
		df<- predict.RNAsnv(cu,df,class="CU")
	}
	#res <- rbind(res.err,res.ai,res.cu)

	df
}


test.models <- function(df,error=NULL,ai=NULL,cu=NULL) {
	if (!is.null(error)) {
		df <- test.RNAsnv(error,df,class="error",response=TRUE)
	}
	if (!is.null(ai)) {
		df<- test.RNAsnv(ai,df,class="AI",response=TRUE)
	}
	if (!is.null(cu)) {
		df<- test.RNAsnv(cu,df,class="CU",response=TRUE)
	}
	#res <- rbind(res.err,res.ai,res.cu)

	df
}



