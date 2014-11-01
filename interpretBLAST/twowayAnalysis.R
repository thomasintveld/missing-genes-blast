
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
# Thomas in't Veld, University of Leuven (last edit 20141101)
# ###########




THRESHOLD <- 45
NROFCORES <- 4
logFile <- 'logFile.log'
LOGFILE <- logFile
library(snow)

load('../data/genelists/genelist.RData')


do.all <- function(){
	zz <- file(logFile, open="wt")
	sink(zz)
	sink(zz, type='message')
	print(paste('beginning analysis', Sys.time()))

	species <- c('gallus', 'anolis', 'taeniopygia', 'latimeria', 'xenopus')
	
	operations <- sapply(species, function(specie) analysisPerSpeciesParallel(specie))

	print(paste('done with primary analysis. Starting interpretation phase.', Sys.time()))

	interpretations <- lapply(species, function(specie) interpretframes(specie))

	print(paste('done with interpretation phase. Tieing all together now.', Sys.time()))

	tieAllTogether()

	print(paste('All done at', Sys.time()))
	sink(type='message')
	sink()

}

analysisPerSpeciesParallel <- function(species){
	
	write(paste('Starting with', species, 'at', Sys.time()), file=LOGFILE, append=FALSE)
	
	filename1 <- paste('../data/results_blast/outputBlast_', species, '_homo.RData', sep='')
	filename2 <- paste('../data/results_blast/outputBlast_homo_', species, '.RData', sep='')

	#filename1 <- paste('homoon',species, '/outputMissingLong.RData',sep='')
	#filename2 <- paste(species, 'onhomo/outputMissingLong.RData', sep='')
	
	# load in these datasets, they both have variable name 'result' so pass them on to different names (df1 and df2)
	load(filename1)
	df1 <- result
	df1$ENSID_REF <- as.character(df1$ENSID_REF)
	df1$ENSID_FOUND <- as.character(df1$ENSID_FOUND)

	load(filename2)	
	df2 <- result
	df2$ENSID_REF <- as.character(df2$ENSID_REF)
	df2$ENSID_FOUND <- as.character(df2$ENSID_FOUND)

	ensids <- unique(df1$ENSID_REF)

	# for testing:
	ensids <- ensids[100:120]

	cl <- makeCluster(NROFCORES)
	clusterExport(cl, c('analysisPerSpecies', 'THRESHOLD', 'interpretframes', 'testAllSplicesForMissing', 'testSpliceForMissing', 'getSplicesFromOneGeneID', 'checkMatchInBothFrames', 'info', 'proteinDictionary', 'log.step', 'LOGFILE' ))
	dfResultTwoWay <- clusterApply(cl, ensids, function(ensid) return(checkMatchInBothFrames(ensid, df1, df2)))
	stopCluster(cl)
	save(dfResultTwoWay, file=paste(species, '_tempResultParallelRun.Rdata', sep=''))
		
	dfResultTwoWay <- do.call('c', dfResultTwoWay)
		
	result <- data.frame(ENSID_REF=ensids, MISSING=dfResultTwoWay)
	save(result, file=paste('../data/results_byprotein/result', species, '.RData', sep=''))
	#write.csv(result, paste(filename1, 'twoway.csv', sep=''), row.names=FALSE)

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

log.step <- function(msg){

	write(paste(msg, 'at', Sys.time()), file=LOGFILE, append=TRUE)	

}

interpretframes <- function(species){

	# get operated list of proteins per species
	filename <- paste('../data/results_byprotein/result', species, '.RData', sep='')
	load(filename)
	df <- result

	# we're using human proteins as reference, so get the links between genes and proteins for homo
	load('../data/genelists/infohomo.RData')
	infohomo <- infosubject

	mergedbig <- merge(df, infohomo, by.x='ENSID_REF', by.y='Ensembl.Protein.ID', all.x=FALSE, all.y=FALSE)

	mergedbig <- unique(mergedbig)

	# Ensembl.Gene.ID is the unique identifier that we have to attain

	ensgids <- unique(mergedbig$Ensembl.Gene.ID) # ~ 20k size, down from 100k

	splicesResult <- testAllSplicesForMissing(mergedbig, ensgids)

	result <- data.frame(Ensembl.Gene.ID=ensgids, Missing=splicesResult)

	save(result, file=paste('../data/results_bygene/', species, '_bygene.RData', sep=''))

	return(result)

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

	species1 <- 'anolis'
	species2 <- 'xenopus'
	species3 <- 'gallus'
	species4 <- 'latimeria'
	species5 <- 'taeniopygia'

	file1 <- paste('../data/results_bygene/', species1, '_bygene.RData', sep='')
	file2 <- paste('../data/results_bygene/', species2, '_bygene.RData', sep='')
	file3 <- paste('../data/results_bygene/', species3, '_bygene.RData', sep='')
	file4 <- paste('../data/results_bygene/', species4, '_bygene.RData', sep='')
	file5 <- paste('../data/results_bygene/', species5, '_bygene.RData', sep='')

	load(file1)
	df1 <- result
	load(file2)
	df2 <- result
	load(file3)
	df3 <- result
	load(file4)
	df4 <- result
	load(file5)
	df5 <- result

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

	result <- df

	fileOut <- paste('../data/results_merged/result_', THRESHOLD, '_', Sys.Date(), '.RData', sep='')

	save(result, file=fileOut)

}

##### HELPERS

proteinDictionary <- function(enspid){

	geneid <- info[info$Ensembl.Protein.ID==enspid,][['Ensembl.Gene.ID']]

	return(geneid)

}

