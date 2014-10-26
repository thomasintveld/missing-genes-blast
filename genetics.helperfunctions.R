blastFind <- function(ensid, fasta, referenceSpecies){

	code <- fasta[[ensid]][[1]]

	log.INFO(paste('Blasting', ensid, 'with ref', referenceSpecies, 'and subject', subjectSpecies, 'at', Sys.time()))

	blastResult <- blastThis(code, referenceSpecies)

	
		result <- interpretBlast(blastResult)
		
		result <- cbind('ENSID_REF'=ensid, result)

	if(WRITEINTERMEDIATETOFILE){

		if(result[1,2]=='MISSING'){
			writeOut(ensid, result[1,2], 'TRUE')
		}else{
			writeOut(ensid, result[1,2], 'FALSE')
		}

	}


	return(result)

}
	
blastThis <- function(seq, referenceSpecies){

	# write sequence
	codeFile <- paste('temp/', as.character(sample(1:10000000,1)), 'code.blast', sep='')
	outFile <- paste('temp/', as.character(sample(1:10000000,1)), 'out.blast', sep='')

	write(seq, file=codeFile)

	# blast against reference species
	# blastformula <- paste("blastp -query ",codeFile," -db ", 
	#	BLASTDBPATH, referenceSpecies, " -out ", outFile, " -outfmt 6 -evalue ", EVALUE, sep="")
	blastformula <- paste(BLASTPATH, "blastp -query ",codeFile," -db ", 
		BLASTDBPATH, referenceSpecies, " -out ", outFile, " -outfmt 6 -evalue ", EVALUE, sep="")
	system(blastformula)

	# # Fields: query id, subject id, % identity, alignment length, 
	#   mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

	blastResult <- try(read.table(outFile, header=FALSE, sep='\t'), silent=TRUE)

	deleteFile(codeFile)
	deleteFile(outFile)

	return(blastResult)


}

deleteFile <- function(file){

	#log.INFO(paste('Deleting file', file))
	try(system(paste("rm", file)), silent=TRUE)

}

interpretBlast <- function(blastResult){

	if(is.null(dim(blastResult))){
		return(MISSINGSTRING)
	}

	blastResult <- blastResult[(blastResult$V3 >= THRESHOLD) & (blastResult$V11 <= EVALUETHRESHOLD),]

	if(dim(blastResult)[[1]] >= 1){

		result <- data.frame('ENSID_FOUND'=blastResult$V2,
							'IDENTITY'=blastResult$V3,
							'EVALUE'=blastResult$V11)

		return(result)
	}else{
		return(MISSINGSTRING)
	}	

}

writeOut <- function(ensid, ensid_ref, result){

	write(paste(ensid, ensid_ref, result, sep=','), file=outputFile, append=TRUE)

}



log.INFO <- function(message){
	write(message, logFile, append=TRUE)

	print(message)

}

loadFasta <- function(file){

	df <- read.fasta(file, as.string=TRUE, forceDNAtolower=FALSE)

	return(df)

}