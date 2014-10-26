
######### STRUCTURE OF ANALYSIS

# 	do.all() : general wrapper function going through the five different species under analysis (using parallel processing)
#	   	|
#	   	|-> 	analysisPerSpecies(specie) 		: analysis function for one species
#	   	|		|_ // load left-to-right missing frame (df1) from previous analysis
#		|		|_ // load right-to-left missing frame (df2) from previous analysis
#		|		|_ checkMatchInBothFrames(ensid, df1, df2) 	: loop over each ENSPID, 
#		|										check if present in other species, check if it links back again. Then it's present, otherwise missing.
#		|
#		|->		interpretframes(specie)			: tie up missing analysis with location data (Chromosome # + Gene.Start.BP) 
#		|												and link to other metadata available (gene name, gc content, length)
#		|
#		[-> 	tieAllTogether()				: do missing algorithm check (currently 'missing' if present in at least three 
#														of (human, anolis, latimeria, xenopus) and missing in both gallus and taeniopygia)
#
# ###########
# Thomas in't Veld, University of Leuven (last edit 20140825)
# ###########




THRESHOLD <- 40
NROFCORES <- 8
logFile <- 'logFile.log'
LOGFILE <- logFile
library(snow)



do.all <- function(){
	zz <- file(logFile, open="wt")
	sink(zz)
	sink(zz, type='message')
	print(paste('beginning analysis on five cores', Sys.time()))

	info <<- bindOneBigInfoFile()

	species <- c('gallus', 'anolis', 'taeniopygia', 'latimeria', 'xenopus')
	
	operations <- sapply(species, function(specie) analysisPerSpeciesParallel(specie))

	print(paste('done with primary analysis. Starting interpretation phase.', Sys.time()))

	interpretations <- clusterApplyLB(cl, species, function(specie) interpretframes(specie))

	print(paste('done with interpretation phase. Tieing all together now.', Sys.time()))

	tieAllTogether()

	print(paste('All done at', Sys.time()))
	sink(type='message')
	sink()

}

analysisPerSpeciesParallel <- function(species){
	
	write(paste('Starting with', species, 'at', Sys.time()), file=LOGFILE, append=FALSE)
	
	filename1 <- paste('homoon',species, '/outputMissingLong.RData',sep='')
	filename2 <- paste(species, 'onhomo/outputMissingLong.RData', sep='')
	
	load(filename1)
	load(filename2)	
	
	ensids <- unique(df1$ENSID_REF)

	# for testing:
	# ensids <- sample(ensids, 200)

	cl <- makeCluster(NROFCORES)
	clusterExport(cl, c('analysisPerSpecies', 'THRESHOLD', 'interpretframes', 'testAllSplicesForMissing', 'testSpliceForMissing', 'getSplicesFromOneGeneID', 'checkMatchInBothFrames', 'info', 'proteinDictionary', 'log.step', 'LOGFILE' ))
	dfResultTwoWay <- clusterApply(cl, ensids, function(ensid) return(checkMatchInBothFrames(ensid, df1, df2)))
	stopCluster(cl)
	save(dfResultTwoWay, file='tempResultParallelRun.Rdata')
		
	dfResultTwoWay <- do.call('c', dfResultTwoWay)
		
	result <- data.frame(ENSID_REF=ensids, MISSING=dfResultTwoWay)
	write.csv(result, paste(filename1, 'twoway.csv', sep=''), row.names=FALSE)

	return(result)

}

save.old.csvs.RData <- function(species){
	
	filename1 <- paste('homoon',species, '/outputMissingLong.csv',sep='')
	filename2 <- paste(species, 'onhomo/outputMissingLong.csv', sep='')
	
	filename1R <- paste('homoon',species, '/outputMissingLong.RData',sep='')
	filename2R <- paste(species, 'onhomo/outputMissingLong.RData', sep='')
	
	df1 <- read.csv(filename1, stringsAsFactors=FALSE)
	df2 <- read.csv(filename2, stringsAsFactors=FALSE)
	
	save(df1, file=filename1R)
	save(df2, file=filename2R)
	

}

analysisPerSpecies <- function(species){
	filename1 <- paste('homoon',species, '/outputMissingLong.csv',sep='')
	filename2 <- paste(species, 'onhomo/outputMissingLong.csv', sep='')

	df1 <- read.csv(filename1, stringsAsFactors=FALSE)
	df2 <- read.csv(filename2, stringsAsFactors=FALSE)

	ensids <- unique(df1$ENSID_REF)

	# for testing:
	# ensids <- sample(ensids, 200)

	dfResultTwoWay <- sapply(ensids, function(ensid) return(checkMatchInBothFrames(ensid, df1, df2)))

	result <- data.frame(ENSID_REF=ensids, MISSING=dfResultTwoWay)
	write.csv(result, paste(filename1, 'twoway.csv', sep=''), row.names=FALSE)

	return(result)

}

checkMatchInBothFrames <- function(ensid, df1, df2){
	
	# log this step
	log.step(ensid)
	
	# 1/ check corresponding ensid in other frame
	matchesLeft <- df1[df1$ENSID_REF==ensid,]
	
		# if nothing there or identity smaller than threshold, return missing
		if(dim(matchesLeft)[[1]]<1){
			return(TRUE)
		}

		# still the case after clearing out non-Threshold identities?
		matchesLeft <- matchesLeft[matchesLeft$IDENTITY > THRESHOLD,]
		if(dim(matchesLeft)[[1]]<1){
			return(TRUE)
		}
	
	# keep only the matchLeft with the largest identity score
	maximumLeft <- max(matchesLeft$IDENTITY)
	matchesLeft <- matchesLeft[matchesLeft$IDENTITY==maximumLeft,]
	matchLeft <- matchesLeft$ENSID_FOUND
			
	# 2/ check if we find the same geneids if we perform the reverse search 
	matchesRight <- df2[df2$ENSID_REF==matchLeft, ]

		# if nothing here or identity smaller than threshold, no linkage and missing!
		if(dim(matchesRight)[[1]]<1){
			return(TRUE)
		}

		# still the case after clearing out non-Threshold identities?
		matchesRight <- matchesRight[matchesRight$IDENTITY > THRESHOLD,]
		if(dim(matchesRight)[[1]]<1){
			return(TRUE)
		}
		
	# keep only the matchRight with the largest identity score
	maximumRight <- max(matchesRight$IDENTITY)
	matchesRight <- matchesRight[matchesRight$IDENTITY==maximumRight,]	
	matchRight <- matchesRight$ENSID_FOUND

	# get all the geneids from initial and found proteinid's (for geneidRight we might have multiple matches with same id)
	geneidLeft <- proteinDictionary(ensid)
	geneidRight <- sapply(matchRight, function(match) return(proteinDictionary(match)))


	# 3/ see if there is at least one mapping between the genes found in Left and the genes found in Right. If so, PRESENT
	if(geneidLeft %in% geneidRight){
		return(FALSE)
	}

	# 4/ ... and else, missing.

	return(TRUE)


}

# checkMatchInThisFrame <- function(ensid, df){

# 	matches <- df[df$ENSID_REF==ensid, ]

# 	if(dim(matches)[[1]]<1){
# 		return(TRUE)
# 	}

# 	if(max(matches$IDENTITY)>THRESHOLD){
# 		return(FALSE)
# 	}else{
# 		return(TRUE)
# 	}

	

# }

# checkMatchInOtherFrame <- function(ensid, df){

# 	matches <- df[df$ENSID_FOUND==ensid, ]

# 	if(dim(matches)[[1]]<1){
# 		return(TRUE)
# 	}

# 	if(max(matches$IDENTITY)>THRESHOLD){
# 		return(FALSE)
# 	}else{
# 		return(TRUE)
# 	}

	

# }


	# Obsolete: old way of calculating twoway analysis
	# # generate missing boolean one way
	# dfResult1 <- sapply(ensids, function(ensid) return(checkMatchInThisFrame(ensid, df1)))
	# # generate missing boolean the other way
	# dfResult2 <- sapply(ensids, function(ensid) return(checkMatchInOtherFrame(ensid, df2)))
	# result <- data.frame(ENSID_REF=ensids, MISSINGLEFT=dfResult1, MISSINGRIGHT=dfResult2)
	# how to verify it's the same gene when present?
	# result['MISSING'] <- ( result$MISSINGLEFT & result$MISSINGRIGHT ) # this is just boolean! need to verify same gene
	

log.step <- function(msg){

	write(paste(msg, 'at', Sys.time()), file=LOGFILE, append=TRUE)	

}

interpretframes <- function(species){
	filename <- paste('homoon',species, '/outputMissingLong.csvtwoway.csv',sep='')
	infohomo <- read.csv('../info/infohomo.csv', stringsAsFactors=FALSE)

	df <- read.csv(filename, stringsAsFactors=FALSE)

	mergedbig <- merge(df, infohomo, by.x='ENSID_REF', by.y='Ensembl.Protein.ID', all.x=FALSE, all.y=FALSE)

	mergedbig <- unique(mergedbig)

	# Ensembl.Gene.ID is the unique identifier that we have to attain

	ensgids <- unique(mergedbig$Ensembl.Gene.ID) # ~ 20k size, down from 100k

	splicesResult <- testAllSplicesForMissing(mergedbig, ensgids)

	result <- data.frame(Ensembl.Gene.ID=ensgids, Missing=splicesResult)

	write.csv(result, paste(species, 'ResultTwoway.csv', sep=''), row.names=FALSE)

	return(result)

}


interpretframes_oneway <- function(species){
	filename <- paste('homoon',species, '/outputMissing.csvtwoway.csv',sep='')
	infohomo <- read.csv('../info/infohomo.csv', stringsAsFactors=FALSE)

	df <- read.csv(filename, stringsAsFactors=FALSE)

	mergedbig <- merge(df, infohomo, by.x='ENSID_REF', by.y='Ensembl.Protein.ID', all.x=FALSE, all.y=FALSE)

	mergedbig <- unique(mergedbig)

	# Ensembl.Gene.ID is the unique identifier that we have to attain

	ensgids <- unique(mergedbig$Ensembl.Gene.ID) # ~ 20k size, down from 100k

	splicesResult <- testAllSplicesForMissing_oneway(mergedbig, ensgids)

	result <- data.frame(Ensembl.Gene.ID=ensgids, Twoway.Missing=splicesResult)

	write.csv(result, paste(species, 'ResultTwoway.csv', sep=''), row.names=FALSE)

	return(TRUE)

}

runner <- function(){

	interpretframes('gallus')
	print('1')
	interpretframes('xenopus')
	print('2')
	interpretframes('latimeria')
	print('3')
	interpretframes('taeniopygia')
	print('4')
	interpretframes('anolis')

}


testAllSplicesForMissing <- function(df, splices){

	result <- sapply(splices, function(splice) testSpliceForMissing(df, splice))

	return(result)

}

testAllSplicesForMissing_oneway <- function(df, splices){

	result <- sapply(splices, function(splice) testSpliceForMissing_oneway(df, splice))

	return(result)

}

testSpliceForMissing <- function(df, gene){

	subset <- getSplicesFromOneGeneID(df, gene)

	#print(paste('Gene', gene, 'spliced at', Sys.time()))

	# missing als alle $MISSING entries TRUE zijn (dan geen enkele match)
	missing <- all(subset$MISSING)

	return(missing)

}

testSpliceForMissing_oneway <- function(df, gene){

	subset <- getSplicesFromOneGeneID(df, gene)

	#print(paste('Gene', gene, 'spliced at', Sys.time()))

	# missing als alle $MISSING entries TRUE zijn (dan geen enkele match)
	missing <- all(subset$MISSING)

	return(missing)

}

getSplicesFromOneGeneID <- function(df, gene){

	# split dataframe to all entries for gene ENSGID
	subset <- df[as.character(df$Ensembl.Gene.ID)==as.character(gene),]

	return(subset)

}

tieAllTogether <- function(){

	file1 <- 'anolisResultTwoway.csv'
	file2 <- 'xenopusResultTwoway.csv'
	file3 <- 'gallusResultTwoway.csv'
	file4 <- 'latimeriaResultTwoway.csv'
	file5 <- 'taeniopygiaResultTwoway.csv'

	df1 <- read.csv(file1, stringsAsFactors=FALSE)
	df2 <- read.csv(file2, stringsAsFactors=FALSE)
	df3 <- read.csv(file3, stringsAsFactors=FALSE)
	df4 <- read.csv(file4, stringsAsFactors=FALSE)
	df5 <- read.csv(file5, stringsAsFactors=FALSE)


	colnames(df1)[2] <- 'anolis'
	colnames(df2)[2] <- 'xenopus'
	colnames(df3)[2] <- 'gallus'
	colnames(df4)[2] <- 'latimeria'
	colnames(df5)[2] <- 'taeniopygia'

	merging1 <- merge(df1, df2, by='Ensembl.Gene.ID')
	merging2 <- merge(df3, df4, by='Ensembl.Gene.ID')
	merging3 <- merge(merging1, df5, by='Ensembl.Gene.ID')
	merging4 <- merge(merging2, merging3, by='Ensembl.Gene.ID')

	df <- merging4

	df['missing1'] <- ((1-df$anolis)*(1-df$xenopus)*(1-df$latimeria)*df$gallus*df$taeniopygia)==1
	df['missing2'] <- ((1-df$anolis)*(1-df$latimeria)*df$gallus*df$taeniopygia)==1
	df['missing3'] <- ((1-df$anolis)*(1-df$xenopus)*df$gallus*df$taeniopygia)==1
	df['missing4'] <- ((1-df$xenopus)*(1-df$latimeria)*df$gallus*df$taeniopygia)==1

	df['missing'] <- (df$missing1 | df$missing2 | df$missing3 | df$missing4)

	fileOut <- paste('resultTwoWayThreshold', THRESHOLD, Sys.time(), '.csv', sep='')

	write.csv(df, fileOut, row.names=FALSE)

}

##### HELPERS

bindOneBigInfoFile <- function(){
	df1 <- normaliseInfoFile(read.csv('../info/infogallus.csv', stringsAsFactors=FALSE))
	df2 <- normaliseInfoFile(read.csv('../info/infoanolis.csv', stringsAsFactors=FALSE))
	df3 <- normaliseInfoFile(read.csv('../info/infolatimeria.csv', stringsAsFactors=FALSE))
	df4 <- normaliseInfoFile(read.csv('../info/infoxenopus.csv', stringsAsFactors=FALSE))
	df5 <- normaliseInfoFile(read.csv('../info/infotaeniopygia.csv', stringsAsFactors=FALSE))
	df6 <- normaliseInfoFile(read.csv('../info/infohomo.csv', stringsAsFactors=FALSE))

	df <- rbind(df1, df2, df3, df4, df5, df6)

	return(df)

}	

normaliseInfoFile <- function(df){

	result <- data.frame(
		Ensembl.Gene.ID=df$Ensembl.Gene.ID,
		Ensembl.Protein.ID=df$Ensembl.Protein.ID,
		Chromosome.Name=df$Chromosome.Name,
		Gene.Start..bp.=df$Gene.Start..bp.,
		Gene.End..bp.=df$Gene.End..bp.,
		Associated.Gene.Name=df$Associated.Gene.Name,
		stringsAsFactors=FALSE
		)
	return(result)
}

proteinDictionary <- function(enspid){

	geneid <- info[info$Ensembl.Protein.ID==enspid,][['Ensembl.Gene.ID']]

	return(geneid)

}

