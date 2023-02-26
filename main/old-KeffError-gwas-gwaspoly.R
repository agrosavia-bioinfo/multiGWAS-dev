#!/usr/bin/Rscript
# INFO   : Script to run GWASpoly for tetraploides (modified from GWASpoly main)
# AUTHOR : Luis Garreta (lgarreta@agrosavia.co)
# DATA   : Feb/2020
# LOG    :
	# r2.0 (Jan/2023): Modified to use GWASpoly 2.3. Removed 2-dom models. Removed plots (no merge plots)
	# r1.1:  Added main. Changed naive kinship matrix to NULL
	# r1.02:  Hidden warnigns qchisq
	# r1.01:  Removed annotations from functions. Output results to "out/" dir 
#-------------------------------------------------------------

main <- function () {
	options (width=300)
	DEBUG  = T
	HOME = Sys.getenv ("MULTIGWAS_HOME")
	library (parallel)
	library (ggplot2)
	library (withr)
	library (tidyr)
	.libPaths (paste0(HOME, "/opt/Rlibs"))
	library (GWASpoly)
	system ("mkdir out")

	args = commandArgs(trailingOnly = TRUE)
	args = c("genotype-gwaspoly-format-example.csv", "phenotype-tuber_shape-example.csv")

	params = list()
	params$genotypeFile      = args [1] 
	params$phenotypeFile     = args [2] 
	params$trait             = colnames (read.csv (params$phenotypeFile))[2]
	params$ploidy            = 4
	params$gwasModel         = "naive"
	#params$snpModels         = c("general","additive","1-dom", "2-dom", "diplo-general", "diplo-additive")
	params$snpModels         = c("1-dom")
	params$correctionMethod  = "Bonferroni"
	params$significanceLevel = 0.05
	params$geneAction        = "dominant"

	msgmsg("Running GWASpoly...")
	a=runToolGwaspoly (params)
 }

set.threshold <- function(data,method="M.eff",level=0.05,n.permute=1000,n.core=1) {
	library (stats)
	library (Matrix)
	library (methods)

	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)
	n.trait <- length(traits)
	models <- colnames(data@scores[[1]])
	model2 <- gsub("-ref","",models,fixed=T)
	model2 <- unique(gsub("-alt","",model2,fixed=T))
	n.model <- length(models)
	methods <- c("M.eff","Bonferroni","FDR","permute")
	stopifnot(is.element(method,methods))
	threshold <- matrix(NA,n.trait,n.model)
	colnames(threshold) <- models
	rownames(threshold) <- traits
	
	if (method=="M.eff") {
	  message ("data@map: ", data@map [,2])
	  chrom <- levels(data@map[,2])
	  message ("chrom: ", chrom)
	  n.chrom <- length(chrom)
	  message ("len chrom: ", n.chrom)
	  r2 <- vector("list",n.chrom)
	  message ("r2: ", r2)
	  names(r2) <- chrom
	  message ("names r2: ", names (r2))
	  for (i in chrom) {
	    ix <- which(data@map[,2]==i)
	    r2[[i]] <- as(cor(data@geno[,ix])^2,"dspMatrix")
	  }
	}
	
	for (i in 1:n.trait) {
		trait <- traits[i]
		if (method=="permute") {
			print(paste("Trait:",trait),quote=F)
			#y <- data@pheno[,trait]
			#ix <- which(!is.na(y))
			max.scores <- matrix(NA,n.permute,n.model)
			colnames(max.scores) <- models
			for (q in 1:n.permute) {
				print(paste("Permutation",q),quote=F)
				data2 <- data
				data2@pheno[,1] <- sample(data@pheno[,1]) #permute id
				data2 <- GWASpoly(data2,models=model2,traits=trait,params=data@params,quiet=T,n.core=n.core)
				for (j in 1:n.model) {max.scores[q,j] <- max(data2@scores[[trait]][,models[j]],na.rm=T)}				
			}
		}
		for (j in 1:n.model) {
			model <-  models[j]
			iv <- which(!is.na(data@scores[[trait]][,model]))
			scores <- as.vector(data@scores[[trait]][iv,model])
			m <- length(scores)
			if (method=="Bonferroni") {threshold[i,j] <- -log10(level/m)}
			if (method=="M.eff") {
			  me <- 0
			  for (chr in chrom) {
			    ix <- data@map[intersect(iv,which(data@map[,2]==chr)),1]
			    message ("IX: ", ix)
			    if (length(ix)>1) {
			      me <- me + Keff(r2=r2[[chr]][ix,ix],alpha=level)
			    } else {
			      me <- me + 1
			    }
			  }
			  threshold[i,j] <- -log10(level/me)
			}
			if (method=="FDR") {
				tmp <- cbind(10^(-scores),.qvalue(10^(-scores)))
				tmp <- tmp[order(tmp[,2]),]
				if (tmp[1,2] > level) {
					threshold[i,j] <- -log10(tmp[1,1])*1.2
				} else {
					k <- max(which(tmp[,2] < level))
					threshold[i,j] <- -log10(mean(tmp[k:(k+1),1]))
				}
			}
			if (method=="permute") {
				threshold[i,j] <- sort(max.scores[,model],decreasing=TRUE)[max(floor(level*n.permute),1)]
			}	
		}
	}
	cat("Thresholds\n")
	print(round(threshold,2))
	return(new("GWASpoly.thresh",map=data@map,pheno=data@pheno,fixed=data@fixed,geno=data@geno,ploidy=data@ploidy,K=data@K,scores=data@scores,effects=data@effects,params=data@params,threshold=threshold))
}

Keff <- function(r2,alpha) {
  message ("Keff r2: ", r2) 
  print (r2) 
  m <- nrow(r2)
  if (m > 1) {
    Q <- sqrt(r2)
    message ("Q: ", Q)
    Q[upper.tri(Q,diag=T)] <- NA
    message ("Keff Q: ", Q) 
    print (class (Q)) 
    save (Q, file="Q.rd") 
    print (class (Q[-1,]))
    print (Q[-1,]) 
    #rmax <- apply(Q[-1,],1,max,na.rm=T)
    rmax <- apply(Q,2,max,na.rm=T)
    print (rmax)
    kappa <- sqrt(1-rmax^(-1.31*log10(alpha)))
    return(1+sum(kappa))
  } else {
    return(1)
  }
}
#-------------------------------------------------------------
# Set parameters and run GWASpoly tool
#-------------------------------------------------------------
runToolGwaspoly <- function (params) {
	NCORES = detectCores ()-1
	# Redirect GAPIT output
	if (DEBUG==F) {
		log <- file("log-GWASpoly-outputs.log", open = "wt")
		sink (log, type="output")
		sink (log, type="message")
	}

	scoresFile = paste0 ("out/tool-GWASpoly-scores-", params$gwasModel, ".csv")

	# Set gene action models (automatic or specific ones)
	# r2.0: Modified by LuisG. Removed 2-dom model as indicated in GWASpoly 2.0
	modelsGwaspoly  = c ("general","additive","diplo-general", "diplo-additive", "1-dom")

	# Default models whether diplo or tetra 
	if (params$geneAction %in% c("all","automatic"))
		snpModes = modelsGwaspoly
	else  # But if user specifided just one model
		if (params$geneAction == "dominant")
			snpModels = c ("1-dom")
		else
			snpModels = c(params$geneAction)

	msgmsg ("SNP models: ",  snpModels)
	params <- append (params, list (snpModels=snpModels))

	# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
	data1 <- read.GWASpoly (pheno.file=params$phenotypeFile, geno.file=params$genotypeFile, 
							ploidy=params$ploidy, format="ACGT", n.traits=1, delim=",")

	# Control population structure
	data2 = controlPopulationStratification (data1, params$gwasModel, NCORES)

	# GWAS execution
	data3 <- runGwaspoly (data2, params, NCORES) 

	# Show results in GWASpoly way
	#msgmsg ("    >>>> GWASpoly showResults...")
	# r2.0: Modified by LuisG. Not mergind manhattan ggplots with default qqplots
	#showResults (data3, params$snpModels, params$trait, params$gwasModel,params$phenotypeFile, params$ploidy)

	# Get SNP associations
	scores  = getQTLGWASpoly (data3, params$gwasModel, params$ploidy)
	colnames (scores)[colnames(scores) %in% c("Chrom","Position")] = c ("CHR","POS")

	scoresColumns = c("MODEL", "GC", "Marker", "CHR", "POS", "P", "SCORE", "THRESHOLD", "DIFF")
	scores <- data.frame (scores[,scoresColumns], scores [,setdiff (colnames(scores), scoresColumns)])

	write.table (file=scoresFile, scores, quote=F, sep="\t", row.names=F)

	# Restore normal output
	if (DEBUG==F) {
		sink ();sink()
	}

	return (list (tool="GWASpoly", scoresFile=scoresFile, scores=scores))
}

#-------------------------------------------------------------
# Control population structure using default Kinship and PCs
#-------------------------------------------------------------
controlPopulationStratification <- function (data1, gwasModel, NCORES) {
	msgmsg ();msgmsg("Controlling populations structure...")

	if (gwasModel=="naive") {
		msgmsg("    >>>> Without any correction") 
		markers       = data1@pheno [,1]
		n             = length (markers)
		kinshipMatrix = matrix (diag (n), n, n, dimnames=list (markers, markers))
		data2         = set.K (data1, LOCO=F, K=NULL, n.core=NCORES)
		gwpParams     = NULL
	}else if (gwasModel == "full") {
		msgmsg("    >>>> Correction using default Kinship and GWASpoly 2.0 LOCO method ")
		kinshipMatrix = NULL
		# r2.0: Modified by LuisG to use the leave-one-chromosome-out (LOCO) method
		N = length (data1@pheno [,1])
		msgmsg ("Population size: ", N)
		data2       = set.K (data1, LOCO=T, n.core=NCORES )
		gwpParams  = set.params (geno.freq=1-5/N, fixed=NULL, fixed.type=NULL)
		write.csv (data2@K,"out/GWASpoly-kinship-matrix.csv", quote=F, row.names=T)
		# r2.0: Modified by LuisG to use maximum genotype frequency threshold instead PCs
		#data2@params  = set.params (n.PC=5, fixed=NULL, fixed.type=NULL)
	}else 
		stop ("Unknown GWAS model (full or naive)")

	return (list (data=data2,gwpParams=gwpParams))
}

#-------------------------------------------------------------
# GWAS execution
#-------------------------------------------------------------
runGwaspoly <- function (data2, params, NCORES) {
	gwasModel        = params$gwasModel
	snpModels        = params$snpModels
	correctionMethod = params$correctionMethod
	signLevel        = params$significanceLevel

 	if (gwasModel %in% c("naive")) {
		msgmsg(">>>> Without params")
		data3 = GWASpoly(data2$data, models=snpModels, traits=NULL, params=NULL, n.core=NCORES)
	}else {
		msgmsg(">>>> With params")
		data3 = GWASpoly(data2$data, models=snpModels, traits=NULL, params=data2$gwpParams, n.core=NCORES)
	}
	
	# QTL Detection
	# r2.0: Modified by LuisG. Set correction method to "M.eff" as default
	data4 = set.threshold (data3, method="M.eff",level=0.05, n.core=NCORES)

	return (data4)
}



#-------------------------------------------------------------
# Plot results
#-------------------------------------------------------------
showResults <- function (data3, models, trait, gwasModel, phenotypeFile, ploidy) {
	msgmsg ();msgmsg("Writing GWASpoly results...")
	#outFile       = paste0 ("out/tool-GWASpoly-scores-", gwasModel)
	#scoresFile    = paste0 (outFile,".csv")
	#scoresFileAll = paste0 (outFile,"-all.csv")
	plotFile      = paste0 ("out/out-GWASpoly-", gwasModel, "-plots.pdf") 

	# QTL Detection
	#data5 = set.threshold (data3, method=correctionMethod,level=0.05,n.core=4)

	# Plot results	
	plotMahattanQQ (plotFile, models, data3, trait, data3, gwasModel, ploidy) 

	#msgmsg(">>>> Writing QTLs to file: ", scoresFile, "...")
	#write.GWASpoly (data5, trait, paste0(scoresFile,".qtls"), "scores", delim="\t")

	#scoresTableAll  = getQTLGWASpoly (data3, gwasModel, ploidy)
	#write.table (file=scoresFileAll, scoresTableAll, quote=F, sep="\t", row.names=F)

	#scoresTableSorted = scoresTableAll [order (scoresTableAll$GC, scoresTableAll$DIFF, decreasing=T),]
	#scoresTable       = scoresTableSorted [!duplicated (scoresTableSorted$Marker, fromLast=F),]
	#write.table (file=scoresFile, scoresTable, quote=F, sep="\t", row.names=F)
}

#-------------------------------------------------------------
# Manhattan and QQ plots
#-------------------------------------------------------------
plotMahattanQQ <- function (plotFile, models, data5, trait, data3, gwasModel, ploidy) {
	msgmsg ("Ploting GWASpoly manhattan...")
	# Create test models for each ref|alt allele if dominant models present
	testModels = c()
	for (m in models)
		if (m=="1-dom") testModels = c (testModels, "1-dom-alt", "1-dom-ref")
		else if (m=="2-dom") testModels = c (testModels, "2-dom-alt", "2-dom-ref")
		else testModels = c (testModels, m)

	n = length (testModels)

	msgmsg ("Plotting models: ", testModels, " into: ", plotFile)

	#pdf (file=plotFile, width=11, height=15)
	#op <- par(mfrow = c(n,2), mar=c(3.5,3.5,2,1), oma=c(0,0,0,0), mgp = c(2.2,1,0)) #MultiGWAS tools
	#for (i in 1:n) {
	#	msgmsg ("Ploting manhattan", i, " ", testModels [i])
	#	manhattan.plot (data3, trait=trait, model=testModels [i])
	#	ggsave (paste0 ("plotM",i,".pdf"))
	#}

	#for (i in 1:n) {
	#	msgmsg ("Ploting QQ ", i, " ", testModels [i])
	#	qqPlot(data3,trait=trait, model=testModels[i], cex=0.3)
	#}

	#plotTitle = sprintf ("%s gwas %s-ploidy for %s trait", gwasModel, ploidy, trait)  
	#mtext(plotTitle, outer=T,  cex=1.5,  line=0)
	#par(op)
	#dev.off()
}

#-------------------------------------------------------------
# Extracts significant QTL
#-------------------------------------------------------------
getQTLGWASpoly <- function(data,gwasModel, ploidy, traits=NULL,models=NULL) {
	stopifnot(inherits(data,"GWASpoly.thresh"))

	if (is.null(traits)) traits <- names(data@scores)
	else stopifnot(is.element(traits,names(data@scores)))

	if (is.null(models)) models <- colnames(data@scores[[1]])
	else stopifnot(is.element(models,colnames(data@scores[[1]])))

	n.model <- length(models)
	n.trait <- length(traits)
	output <- data.frame(NULL)
	for (j in 1:n.model) {
		#ix <- which(data@scores[[traits[1]]][,models[j]] > (data@threshold[traits[1],models[j]]) - 1)
		ix <- which (data@scores[[traits[1]]][,models[j]] != 0)
		markers <-  data.frame (SNP=data@map[ix,c("Marker")])

		scores <- data@scores[[1]][,models[j]]
		datax = calculateInflationFactor (scores)

		n.ix <- length(ix)
		
		gc=rep(datax$delta,n.ix) 
		scores=round(data@scores[[traits[1]]][ix,models[j]],2)
		thresholds=round(rep(data@threshold[traits[1],models[j]],n.ix),2)
		diffs = (scores - thresholds)
		pvalues = 10^(-scores)
		df = data.frame(Ploidy=rep (ploidy, n.ix), Type=rep (gwasModel, n.ix),
						data@map[ix,], GC=gc, MODEL=rep(models[j],n.ix),
						P=pvalues,SCORE=scores, THRESHOLD=thresholds, DIFF=diffs,
						Effect=round(data@effects[[traits[1]]][ix,models[j]],2))
						#stringsAsFactors=F,check.names=F)

		output <- rbind(output, df)
	}
	#out <-cbind (Type=gwasModel, output)
	#output <- output [order(-output$GC, -output$DIFF),]
	output = output [order(-output$DIFF),]
	#output = output [!duplicated (output$Marker),]
	#outputPositives = output [output$DIFF > 0,]
	#outputNegatives = output [output$DIFF <= 0,]

	#outQTLsAllSNPs = rbind (outputPositives, outputNegatives)

	return(output)
}

#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
# It can fire warning, here they are hidign
#-------------------------------------------------------------
calculateInflationFactor <- function (scores) {
	oldw <- getOption("warn")
	options(warn = -1)

	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	options (warn = oldw)

	return (list(delta=delta, scores=x))
}

#-------------------------------------------------------------
# QQ plot
#-------------------------------------------------------------
qqPlot <- function(data,trait,model,cex=1,filename=NULL) {
	stopifnot(inherits(data,"GWASpoly.fitted"))
	traits <- names(data@scores)

	stopifnot(is.element(trait,traits))
	models <- colnames(data@scores[[trait]])
	stopifnot(is.element(model,models))
	scores <- data@scores[[trait]][,model]

	datax = calculateInflationFactor (scores)

	n <- length(datax$scores)
	unif.p <- -log10(ppoints(n))
	if (!is.null(filename)) {postscript(file=filename,horizontal=FALSE)}
	par(pty="s")
	plot(unif.p, datax$scores, pch=16,cex=cex,
		 xlab=expression(paste("Expected -log"[10],"(p)",sep="")),
		 ylab=expression(paste("Observed -log"[10],"(p)",sep="")),
		 main=paste(trait," (",model,") ",sep=""))

	mtext (bquote(lambda[GC] == .(datax$delta)), side=3, line=-2, cex=0.7)

	lines(c(0,max(unif.p)),c(0,max(unif.p)),lty=2)
	if (!is.null(filename)) {dev.off()}
	return(datax$delta)
}

#-------------------------------------------------------------
# Read input genotype and genotype (format: "numeric", "AB", or "ACGT")
#-------------------------------------------------------------
initGWAS <- function (phenotypeFile, genotypeFile, ploidy, format="ACGT") {
	msgmsg ();msgmsg("Initializing GWAS...");msgmsg ()

	data1 <- read.GWASpoly (ploidy = ploidy, pheno.file = phenotypeFile, 
							geno.file = genotypeFile, format = "ACGT", n.traits = 1, delim=",")

	return (data1)
}
#-------------------------------------------------------------
# Add label to filename
#-------------------------------------------------------------
addLabel <- function (filename, label)  {
	nameext = strsplit (filename, split="[.]")
	newName = paste0 (nameext [[1]][1], "-", label, ".", nameext [[1]][2])
	return (newName)
}
#-------------------------------------------------------------
#-------------------------------------------------------------
msgmsg <- function (...) {
  messages = unlist (list (...))
  cat ("\t>>>>", messages, "\n")
}

#-------------------------------------------------------------
#main ()
