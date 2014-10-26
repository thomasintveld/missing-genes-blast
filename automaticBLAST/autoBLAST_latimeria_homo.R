require('seqinr')
require('snow')

BLASTPATH <- "blast/bin/"
BLASTDBPATH <- "blast/db/"
FASTAHOMO <- "species/Homo_sapiens.GRCh37.75.pep.all.fa"
FASTAGALLUS <- "species/Gallus_gallus.Galgal4.75.pep.all.fa"
FASTATAENIOPYGIA <- "species/Taeniopygia_guttata.taeGut3.2.4.75.pep.all.fa"
FASTALATIMERIA <- "species/Latimeria_chalumnae.LatCha1.75.pep.all.fa"
FASTAXENOPUS <- "species/Xenopus_tropicalis.JGI_4.2.75.pep.all.fa"

#### PUT SPECIES TO ANALYSE HERE

referenceSpecies <- "latimeria"
subjectSpecies <- "homo"
fastaFile <- FASTAHOMO

####

# referenceSpecies<- "gallus"
# referenceSpecies <- "taeniopygia"
# referenceSpecies <- "anolis"
# referenceSpecies <- "xenopus"
# referenceSpecies <- "latimeria"


THRESHOLD <- 45
EVALUE <- 1E-5
EVALUETHRESHOLD <- 0.0001
MISSINGSTRING <- data.frame("ENSID_FOUND"='MISSING', "IDENTITY"=0, "EVALUE"=0)
NROFCORES <- 8
PARALLEL <- TRUE
TESTRUN <- FALSE
WRITEINTERMEDIATETOFILE <- FALSE


outputFile <- paste('outputMissing','_', referenceSpecies, '_', subjectSpecies, '_',Sys.Date(), '.csv', sep='')
resultFileRData <- paste('outputBlast','_', referenceSpecies, '_', subjectSpecies, '_',Sys.Date(), '.RData', sep='')
codeFile <- paste('code','_', referenceSpecies, '_', subjectSpecies, sep='')
outFile <- paste('out','_', referenceSpecies, '_', subjectSpecies, sep='')
logFile <- paste('log','_', referenceSpecies, '_', subjectSpecies, '_',Sys.Date(),'.log', sep='')

source('genetics.helperfunctions.R')

# The idea is to blast genes from a 'subject' species with a 'reference' species.
# We load in the genetic information for the subject species with both a CSV and a FASTA file;
# the CSV contains a list of the needed proteins; the FASTA contains the protein sequence info.

if(!exists('infosubject')){
	infosubject <- read.csv(paste('info/info', subjectSpecies, '.csv', sep=''), 
								header=TRUE, stringsAsFactors=FALSE)
}

speciesENS <- as.character(infosubject$Ensembl.Protein.ID)

# Load FASTA for subject

if(!exists('fastaFromFile')){
	log.INFO("Loading FASTA for subject species, this takes a bit of time.")
	fastaFromFile <- loadFasta(fastaFile)
}

fasta <- fastaFromFile

# If we are in a test-run, load small subset of genetic data

if(TESTRUN){
	allensid <- sample(speciesENS, 100)
}else{
	allensid <- speciesENS
}

if(PARALLEL){
	cl <- makeCluster(NROFCORES)
	clusterExport(cl, list=ls())

}

# Run this as the main function, e.g. df <- main()

main <- function(){
	
	if(WRITEINTERMEDIATETOFILE){
		# Initalise the output file
		write(paste('ENSID_REF', 'ENSID_FOUND', 'MISSING', sep=","), file=outputFile, append=FALSE)
	}

	write('Starting all', logFile, append=FALSE)
	log.INFO(paste('started processing at', Sys.time()))

	if(PARALLEL){
		clusterExport(cl, list=ls())
		clusterExport(cl, c('blastFind', 'blastThis', 'interpretBlast'))
		result <- clusterApplyLB(cl, allensid, function(ensid) blastFind(ensid, fasta, referenceSpecies))
		stopCluster(cl)
	}else{
		result <- lapply(allensid, function(ensid) blastFind(ensid, fasta, referenceSpecies))

	}

	result <- do.call('rbind', result)

	log.INFO(paste('Ended all processing at', Sys.time()))
	
	# save result to RData
 	save(result, file=resultFileRData)
 	
 	return(result)
}




