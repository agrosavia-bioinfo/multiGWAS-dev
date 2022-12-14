#!/usr/bin/Rscript

source ("lglib10.R")
library(dplyr)

#-----------------------------------------------------------
# Sort best SNPs from multiple tools using a own score
#-----------------------------------------------------------
sortBestSNPsbyMGScore <- function (scoresFilename, nBest, tool, geneAction) {
	scoresTable = read.table (scoresFilename, sep="\t", header=T)
	view (scoresTable, n=0,m=0)

	# Add Count of SNPs between groups
	dfCountSNPs  = data.frame (add_count (scoresTable, SNP, sort=F, name="scoreShared")); 

	# Score GC
	scoreSign   = ifelse (scoresTable$SIGNIFICANCE, 1,0)
	scoreGC     = 1 - abs (1-scoresTable$GC)
	scoreSNPs   = data.frame (add_count (scoresTable, SNP))$n

	SCORE_SNP = 0.5*scoreGC + 0.3*scoreSNPs + 0.2*scoreSign

	dataSNPsScores  = select (data.frame (scoresTable, SCORE_SNP), SNP, SCORE_SNP, everything())
	dataSNPsScores  = dataSNPsScores [order (-SCORE_SNP),] 
	write.csv (dataSNPsScores, "out-BestSNPs-ScoreMG.csv", row.names=F)

