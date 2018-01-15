# mutantAligner.R - Last updated: 01-09-18
# Developed By: Bob Morrison
# Sequence deconvolution code for Transcriptional Regulon Induced Phenotype (TRIP) screening. 
# Generates relative abundance information from raw fastq sequencing files
# This code depends on the DuffyNGS and DuffyTools packages, which can be found at:
# https://github.com/sturkarslan/DuffyNGS and https://github.com/sturkarslan/DuffyTools

library( DuffyNGS)
require( Biostrings)
multicore.setup(4)
setCurrentSpecies( "MT_H37")

# constants to folders & files
BOWTIE2_PROGRAM <- "./bin/bowtie2"

TFOE.BowtieIndex <- "./BowtieIndexes/MTb.TFOE.MutantGenome"
TFOE.GenomeFile <- "./CustomMutantGenomes/MTb.TFOE.MutantGenome.fasta"
TFOE.CompleteGeneFile <- "./CustomMutantGenomes/TFOE.Complete207.Primers.txt"
TFKO.BowtieIndex <- "./BowtieIndexes/MTb.TFKO.MutantGenome"
TFKO.GenomeFile <- "./CustomMutantGenomes/MTb.TFKO.MutantGenome.fasta"
TFKO.CompleteGeneFile <- "./CustomMutantGenomes/TFKO.Primers.txt"

SampleKeyFolder <- "SampleKeys"
ResultsFolder <- "Results"
SampleKeyColumns <- c("SampleID", "Experiment", "Day", "Replicate", "Induced", "Class",
				"SubClass1", "SubClass2", "MutantPool", "File1", "File2", "FastqPath")
MetaDataColumns <- SampleKeyColumns[ 1:9]
ModelMetaDataColumns <- c( MetaDataColumns, "Run", "ReadCount")

DAY_ZERO <- c( 0, 1)
MIN_LOG2_RPMHEG <- 4

checkX11( bg="white", width=10, height=7)
par( mai=c(0.6,0.4,0.6, 0.4))


# function for the standard TFOE Transcription Factor Over Expression system runs
# align, count up pairs, assess RPM, and write out count data for all samples
align.OneSample <- function( sampleKeyFile="SamplesKey.Run1.47TFlogPre.txt", doBowtie=FALSE, dropZeroGenes=TRUE, 
				bowtieMode=c("Bowtie2","Bowtie1"), verbose=FALSE, doPairs=FALSE, makePies=TRUE, 
				mode=c("TFOE","TFKO"), makeROC=FALSE, doMODEL=TRUE, sharedDayZeroUninduced=NULL,
				trim5=0, trim3=0) {
	
	# set the needed Genome & Target index for TFOE
	mode <- match.arg( mode)
	if ( mode == "TFOE") {
		BowtieIndex <<- TFOE.BowtieIndex
		GenomeFile <<- TFOE.GenomeFile
		CompleteGeneFile <<- TFOE.CompleteGeneFile
	}
	if ( mode == "TFKO") {
		BowtieIndex <<- TFKO.BowtieIndex
		GenomeFile <<- TFKO.GenomeFile
		CompleteGeneFile <<- TFKO.CompleteGeneFile
	}

	# get all the Rv Gene names
	MTgenome <- loadFasta( GenomeFile, verbose=F)
	genes <- grep( "^Rv|Empty", MTgenome$desc, value=T)
	nGenes <- length(genes)
	cat( "\nN_Genes in MTb Mutant Genome file:    ", nGenes)

	# determine which alignment tool to use
	bowtieMode <- match.arg( bowtieMode) 
	cat( "\nUsing program:   ", toupper( bowtieMode), "\n")

	# load that sample key
	allSamples <- loadSampleKeyFile( sampleKeyFile)
	nSamples <- nrow( allSamples)

	# the 'Run Name' is extracted from the name of the Samples Key filename for the results folder name
	runName <- extractRunName( sampleKeyFile)
	results.path <- createResultsFolder( sampleKeyFile)

	# some storage is for the entire run... even when more than one experiment is in the Sample Key file
	outStats <<- matrix( nrow=nSamples, ncol=8)
	rownames(outStats) <<- allSamples$SampleID
	colnames(outStats) <<- c( "StartingReads", "ValidMate1", "Pct_ValidMate1", "ValidMate2", 
				"Pct_ValidMate2", "ValidPairs", "Pct_ValidPairs", "N_Genes")

	# let's allow more than one 'Experiment' per Sample Key file
	expFactor <- factor( allSamples$Experiment)
	tapply( 1:nrow(allSamples), expFactor, function(x) {

		# do the subset of samples for this experiment
		samples <- allSamples[x, ]
		nSamples <- nrow( samples)

		# get the details about this experiment
		myExperiment <- samples$Experiment[1]
		cat( "\n\nProcessing Experiment:   ", myExperiment)

		# make a big matrix to hold all the data
		out <- matrix( 0, nrow=nGenes, ncol=nSamples)
		colnames(out) <- samples$SampleID
		rownames(out) <- genes
		poolNames <- nExpectedGenes <- vector( length=nSamples)
		#poolLists <- vector( mode="list")

		# let's do the alignments on multiple cores at once
		doOneMutantAlign <- function( i) {
			sid <- samples$SampleID[i]
			doMutantAlignments( sid, samples$File1[i], samples$File2[i], fastqPath=samples$FastqPath[i],
					doBowtie=doBowtie, bowtieMode=bowtieMode, results.path=results.path, 
					trim5=trim5, trim3=trim3, verbose=verbose)
		}
		multicore.lapply( 1:nSamples, FUN=doOneMutantAlign)
		
		# next, let's tabulate all the read pairs, for each sample on multiple cores at once...
		TABULATE.FUN <- tabulateMutantAlignments.TFOE
		if (mode == "TFKO") TABULATE.FUN <- tabulateMutantAlignments.TFKO
		multicore.lapply( samples$SampleID, FUN=TABULATE.FUN, doBowtie=doBowtie, bowtieMode=bowtieMode,
				doPairs=doPairs, results.path=results.path)

		# now we can do the final summaary and pie plots
		for ( i in 1:nSamples) {
			sid <- samples$SampleID[i]
			# stick the alignment results into the overall run table, not just for this one experiment
			rowptr <- x[i]
			ans <- NULL
			if ( mode == "TFOE") {
				ans <- summarizeMutantAlignments.TFOE( sid, bowtieMode=bowtieMode, rowptr=rowptr, 
					results.path=results.path, makePie=makePies)
			} else if (mode == "TFKO") {
				ans <- summarizeMutantAlignments.TFKO( sid, bowtieMode=bowtieMode, rowptr=rowptr, 
					results.path=results.path, makePie=makePies)
			}
			if ( is.null(ans)) next

			# put these counts where they go for this experiment
			where <- match( names(ans), genes, nomatch=0)
			out[ where, i] <- ans[ where > 0]

			# get the details about this pool
			myPool <- samples$MutantPool[i]
			poolGenes <- getMutantPoolGenes( myPool)
			nGenesThisPool <- length( poolGenes)
			poolNames[i] <- myPool
			nExpectedGenes[i] <- length(poolGenes)
		}

		if (makeROC) {
			ROC.rawCount.cutoff <- makeROCimage( out, poolGenes=poolGenes, text="Raw Counts", 
							experiment=myExperiment, bowtieMode=bowtieMode, 
							results.path=results.path)
		} else {
			ROC.rawCount.cutoff <- 10
		}

		sumPerGene <- apply( out, MARGIN=1, FUN=sum)
		nDetect <- sum( sumPerGene > ROC.rawCount.cutoff)
		cat( "\n\nN_Detected Genes over", nSamples, "samples:  ", nDetect, "\tROC cut: ", ROC.rawCount.cutoff)
		if ( nDetect < nGenes && dropZeroGenes) {
			toDrop <- which( sumPerGene == 0)
			out <- out[ -toDrop, ]
			cat( "\nDropped genes with zero reads:  ", length( toDrop))
		}

		# normalize to RPM and write the various result files
		writeResultTables( out, samples, results.path=results.path, nExpectedGenes=nExpectedGenes, 
					poolGenes=poolGenes, experiment=myExperiment, bowtieMode=bowtieMode,
					makeROC=makeROC, runName=runName) 

		cat( "\nFinished Experiment:  ", myExperiment, "\n")

	})  ## done with tapply() of all experiments in the run

	# write out the alignment stats overview table too..
	os <- data.frame( "SampleID"=rownames(outStats), outStats, stringsAsFactors=F)
	outfile <- file.path( results.path, paste( "AlignmentStatistics", bowtieMode, "txt", sep="."))
	write.table( os, outfile, sep="\t", quote=F, row.names=F)

	plotAlignmentOverview( outStats, allSamples, bowtieMode, results.path=results.path)

	# when we are all done, try to run the growth defect code as a last step
	if ( doMODEL && mode == "TFOE") model.OneTFOEsample( sampleKeyFile, sharedDayZeroUninduced=sharedDayZeroUninduced)
	if ( doMODEL && mode == "TFKO") model.OneTFKOsample( sampleKeyFile)
	
	cat( "\nDone.\n")
}


# function to do Bowtie alignment for both mate files of one sample
doMutantAlignments <- function( sid, file1, file2, fastqPath=".", doBowtie=TRUE, 
				bowtieMode=c("Bowtie2", "Bowtie1"), results.path=".", 
				trim5=0, trim3=0, verbose=TRUE) {

	# verify FASTQ files are found
	file1 <- file.path( fastqPath, file1)
	file2 <- file.path( fastqPath, file2)
	if ( ! file.exists( file1)) {
		cat("\n\nError:  Mate 1 Fastq file not found: ", sid, "\tFile: ", file1)
		stop()
	}
	if ( ! file.exists( file2)) {
		cat("\n\nError:  Mate 2 Fastq file not found: ", sid, "\tFile: ", file2)
		stop()
	}

	bowtieMode <- match.arg( bowtieMode)
	bamPath <- file.path( results.path, paste( "BAMresults", bowtieMode, sep="."))
	if ( ! file.exists( bamPath)) dir.create( bamPath, recursive=T)
	outfile1 <- file.path( bamPath, paste( sid, 1, "bam", sep="."))
	outfile2 <- file.path( bamPath, paste( sid, 2, "bam", sep="."))

	if ( ! doBowtie && all( file.exists( c( outfile1, outfile2)))) {
		cat( "\nSkip..  BAM files already exist for: ", sid)
		return()
	}

	cat( "\n\n\n-----------------------------------------------------------------\n\nCalling BOWTIE alignment on mutant sample:  ", sid)

	# build the Unix command line
	trimOptions <- ""
	if (trim5 > 0) trimOptions <- paste( trimOptions, " --trim5", trim5)
	if (trim3 > 0) trimOptions <- paste( trimOptions, " --trim3", trim3)

	if ( bowtieMode == "Bowtie1") {
		program <- Sys.getenv( "BOWTIE_PROGRAM")
		options <- paste( " -q -B 1 -v 2 --best -M 1 --sam --threads 4 --ignore-quals ", trimOptions)
		fileFlag <- ""
	} else {
		program <- Sys.getenv( "BOWTIE2_PROGRAM")
		program <- BOWTIE2_PROGRAM
		options <- paste( " -q  --local  --score-min G,20,9  --very-sensitive-local  --threads 4 --ignore-quals ",
					trimOptions, " -x ")
		fileFlag <- " -U "
	}

	# catch the statistics so we can echo them
	catchStdErrFile <- paste( sid, "bowtieStdErr.txt", sep=".")
	catchStdErrLog <- paste( "  2> ", catchStdErrFile)

	cmdline <- paste( program, options, BowtieIndex, fileFlag, file1, catchStdErrLog, " | samtools view -b - >", outfile1)
	cat( "\n\nDoing Mate 1: ", sid, "\n")
	if (verbose) cat( "\nBowtie command line: \n", cmdline, "\n")
	file.delete( catchStdErrFile)
	system( cmdline)
	lines <- readLines( catchStdErrFile)
	writeLines( lines)

	cmdline <- paste( program, options, BowtieIndex, fileFlag, file2, catchStdErrLog, " | samtools view -b - >", outfile2)
	cat( "\nDoing Mate 2: ", sid, "\n")
	if (verbose) cat(  "\nBowtie command line: \n", cmdline, "\n")
	file.delete( catchStdErrFile)
	system( cmdline)
	lines <- readLines( catchStdErrFile)
	writeLines( lines)

	file.delete( catchStdErrFile)
	if (verbose) cat( "\nAlignment Done:  ", sid, "\n")
	return()
}


# function to scan the 2 BAM files of TFOE data and tabulate the results by combining the mate pairs
tabulateMutantAlignments.TFOE <- function( sid, bowtieMode="Bowtie2", results.path=".", 
				makePie=TRUE, doBowtie=TRUE, doPairs=TRUE, trackNoHits=TRUE) {

	N_NOHIT_KEEP <- 100

	# where the BAM files will be found
	bamPath <- file.path( results.path, paste( "BAMresults", bowtieMode, sep="."))
	infile1 <- paste( sid, 1, "bam", sep=".")
	infile2 <- paste( sid, 2, "bam", sep=".")
	infile1 <- file.path( bamPath, infile1)
	infile2 <- file.path( bamPath, infile2)

	# where the results will be put
	pairPath <- file.path( results.path, paste( "PAIRresults", bowtieMode, sep="."))
	if ( ! file.exists( pairPath)) dir.create( pairPath, recursive=T)
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))

	# let's try to find/call empty adapters too
	adapter1 <- "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapter2 <- "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	rcadapter1 <- myReverseComplement( adapter1)
	rcadapter2 <- myReverseComplement( adapter2)
	adapterMinScore <- 30

	# do the tabulation whenever Bowtie was run or the answer not already made
	doTabulate <- ( doBowtie || doPairs || !file.exists( pairFile))
	if ( doTabulate) {

		cat( "\n\nReading BAM alignment results for sample:  ", sid, "  ")
		# set up to iterate thru the two files
		reader1 <- bamReader( infile1)
		reader2 <- bamReader( infile2)
		refData <- getRefData( reader1)
	
		hasMore <- TRUE
		big1 <- big2 <- bigID <- vector()
		big1pos <- big2pos <- vector()
		nohitTbl1 <- nohitTbl2 <- NULL

		while (hasMore) {
			chunk1 <- getNextChunk( reader1, n=100000)
			chunk2 <- getNextChunk( reader2, n=100000)
			nNow <- size( chunk1)
			if ( nNow < 100000) hasMore <- FALSE
			if ( nNow < 1) break
	
			readIDs1 <- readID( chunk1)
			readIDs2 <- readID( chunk2)
			# some Sequencing cores will leave mate suffixes...  strip them off
			readIDs1 <- sub( "/1$", "", readIDs1)
			readIDs2 <- sub( "/2$", "", readIDs2)
			refIDs1 <- refID( chunk1)
			refIDs2 <- refID( chunk2)
			seqIDs1 <- refID2seqID( refIDs1, refData=refData)
			seqIDs2 <- refID2seqID( refIDs2, refData=refData)
			seqIDs1[ is.na( seqIDs1)] <- "NoHit"
			seqIDs2[ is.na( seqIDs2)] <- "NoHit"
			pos1 <- position( chunk1)
			pos2 <- position( chunk2)
			isNoHit1 <- which( seqIDs1 == "NoHit")
			if ( length( isNoHit1)) {
				seqDNA1 <- readSeq(chunk1)
				# check for empty adaptor ends
				scores <- pairwiseAlignment( seqDNA1[isNoHit1], rcadapter2, type='local', scoreOnly=TRUE)
				reLabel <- isNoHit1[ scores >= adapterMinScore]
				if (length(reLabel)) seqIDs1[ reLabel] <- "EmptyAdapter"
				if (trackNoHits) {
					smlTbl <- sort( table( seqDNA1[ isNoHit1]), decreasing=T)
					Nkeep <- min( N_NOHIT_KEEP, length(smlTbl))
					nohitTbl1 <- if (is.null(nohitTbl1)) smlTbl[ 1:Nkeep] else mergeTables( nohitTbl1, smlTbl[1:Nkeep])
				}
				# check for Poly N
				isNoHit1 <- setdiff( isNoHit1, reLabel)
				bases <- strsplit( seqDNA1[isNoHit1], split="")
				nPolyN <- sapply( bases, function(x) sum( x == "N"))
				isPolyN <- which( nPolyN > 5)
				if (length(isPolyN)) seqIDs1[ isNoHit1[isPolyN]] <- "PolyN"
			}
			isNoHit2 <- which( seqIDs2 == "NoHit")
			if ( length( isNoHit2)) {
				seqDNA2 <- readSeq(chunk2)
				# check for empty adaptor ends
				scores <- pairwiseAlignment( seqDNA2[isNoHit2], rcadapter1, type='local', scoreOnly=TRUE)
				reLabel <- isNoHit2[ scores >= adapterMinScore]
				if (length(reLabel)) seqIDs2[ reLabel] <- "EmptyAdapter"
				if (trackNoHits) {
					smlTbl <- sort( table( seqDNA2[ isNoHit2]), decreasing=T)
					Nkeep <- min( N_NOHIT_KEEP, length(smlTbl))
					nohitTbl2 <- if (is.null(nohitTbl2)) smlTbl[ 1:Nkeep] else mergeTables( nohitTbl2, smlTbl[1:Nkeep])
				}
				# check for Poly N
				isNoHit2 <- setdiff( isNoHit2, reLabel)
				bases <- strsplit( seqDNA2[isNoHit2], split="")
				nPolyN <- sapply( bases, function(x) sum( x == "N"))
				isPolyN <- which( nPolyN > 5)
				if (length(isPolyN)) seqIDs2[ isNoHit2[isPolyN]] <- "PolyN"
			}
	
			# force the mate pairs to match by their IDs
			if ( ! all( readIDs1 == readIDs2)) {
				both <- intersect( readIDs1, readIDs2)
				wh1 <- match( both, readIDs1)
				wh2 <- match( both, readIDs2)
				refIDs1 <- refIDs1[wh1]
				refIDs2 <- refIDs2[wh2]
				seqIDs1 <- seqIDs1[wh1]
				seqIDs2 <- seqIDs2[wh2]
				pos1 <- pos1[wh1]
				pos2 <- pos2[wh2]
				readIDs1 <- readIDs2 <- both
			}
			# sanity check that we have intact pairs
			if ( length( readIDs1) < (nNow/2)) stop( "Error:  Problem pairing up Mate Read IDs.")
			cat( ".")
			big1 <- c( big1, seqIDs1)
			big2 <- c( big2, seqIDs2)
			bigID <- c( bigID, readIDs1)
			big1pos <- c( big1pos, pos1)
			big2pos <- c( big2pos, pos2)
		}
	
		bamClose( reader1)
		bamClose( reader2)
		cat( "Done.\n")

		cat( "\nForce mate pairs into alphabetical order..")
		# make the pairs be in alpha order
		# we are using the fact that "Anchor" comes before all the "Rv" genes when we sort by alphabetical order
		# "Anchor" will be in 'big1' and the Rv gene will be in 'big2'
		Npairs <- length(big1)
		isAlpha <- (big1 < big2)
		isOne <- rep.int( 1, Npairs)
		isTwo <- rep.int( 2, Npairs)
		if ( any( ! isAlpha)) {
			toSwap <- which( ! isAlpha)
			isOne[toSwap] <- 2
			isTwo[toSwap] <- 1
		}
			
		# each entry should be one Anchor and one Gene
		cat( "\nAssess pairing success..")
		text1 <- text2 <- rep.int( "GENE", Npairs)
		text1[ grep( "Anchor", big1, fixed=T)] <- "ANCHOR"
		text2[ grep( "Anchor", big2, fixed=T)] <- "ANCHOR"
		text1[ grep( "NoHit", big1, fixed=T)] <- "NOHIT"
		text2[ grep( "NoHit", big2, fixed=T)] <- "NOHIT"
		text1[ grep( "EmptyAdapter", big1, fixed=T)] <- "EMPTY_ADAPT"
		text2[ grep( "EmptyAdapter", big2, fixed=T)] <- "EMPTY_ADAPT"
		text1[ grep( "PolyN", big1, fixed=T)] <- "POLY_N"
		text2[ grep( "PolyN", big2, fixed=T)] <- "POLY_N"
		bigStr <- paste( text1, text2, sep=".")

		pairDF <- data.frame( "ReadID"=bigID, "Mate1"=big1, "Mate2"=big2, "Status"=bigStr, 
				"Pos1"=big1pos, "Pos2"=big2pos, stringsAsFactors=FALSE)
		save( pairDF, file=pairFile)
		cat( "\nWrote Read Pair results:  ", sid, "\tN_Pairs: ", Npairs)

		if (trackNoHits && !is.null( nohitTbl1)) {
			nohitFile <- file.path( pairPath, paste( sid, "TopNoHits.Mate1.fasta", sep="."))
			nohitTbl1 <- sort( nohitTbl1, decreasing=T)
			if ( length(nohitTbl1) > N_NOHIT_KEEP) nohitTbl1 <- nohitTbl1[ 1:N_NOHIT_KEEP]
			NFA <- length(nohitTbl1)
			writeFasta( as.Fasta( desc=paste( "Mate1_", 1:NFA, "_Count=", nohitTbl1, sep=""),
					seq=names(nohitTbl1)), nohitFile, line.width=200)
		}
		if (trackNoHits && !is.null( nohitTbl2)) {
			nohitFile <- file.path( pairPath, paste( sid, "TopNoHits.Mate2.fasta", sep="."))
			nohitTbl2 <- sort( nohitTbl2, decreasing=T)
			if ( length(nohitTbl2) > N_NOHIT_KEEP) nohitTbl2 <- nohitTbl2[ 1:N_NOHIT_KEEP]
			NFA <- length(nohitTbl2)
			writeFasta( as.Fasta( desc=paste( "Mate2_", 1:NFA, "_Count=", nohitTbl2, sep=""),
					seq=names(nohitTbl2)), nohitFile, line.width=200)
		}

	} else {
		cat( "\nSkip..  Read Pair results exist already: ", sid)
	}
}


# function to scan the 2 BAM files of TFKO data and tabulate the results by combining the mate pairs
tabulateMutantAlignments.TFKO <- function( sid, bowtieMode="Bowtie2", results.path=".", 
				makePie=TRUE, doBowtie=TRUE, doPairs=TRUE, trackNoHits=TRUE) {

	N_NOHIT_KEEP <- 100

	# where the BAM files will be found
	bamPath <- file.path( results.path, paste( "BAMresults", bowtieMode, sep="."))
	infile1 <- paste( sid, 1, "bam", sep=".")
	infile2 <- paste( sid, 2, "bam", sep=".")
	infile1 <- file.path( bamPath, infile1)
	infile2 <- file.path( bamPath, infile2)

	# where the results will be put
	pairPath <- file.path( results.path, paste( "PAIRresults", bowtieMode, sep="."))
	if ( ! file.exists( pairPath)) dir.create( pairPath, recursive=T)
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))

	# let's try to find/call empty adapters too
	adapter1 <- "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	adapter2 <- "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	rcadapter1 <- myReverseComplement( adapter1)
	rcadapter2 <- myReverseComplement( adapter2)
	adapterMinScore <- 30

	# do the tabulation whenever Bowtie was run or the answer not already made
	doTabulate <- ( doBowtie || doPairs || !file.exists( pairFile))
	if ( doTabulate) {

		cat( "\n\nReading BAM alignment results for sample:  ", sid, "  ")
		# set up to iterate thru the two files
		reader1 <- bamReader( infile1)
		reader2 <- bamReader( infile2)
		refData <- getRefData( reader1)
	
		hasMore <- TRUE
		big1 <- big2 <- bigID <- vector()
		big1pos <- big2pos <- vector()
		nohitTbl1 <- nohitTbl2 <- NULL

		while (hasMore) {
			chunk1 <- getNextChunk( reader1, n=100000)
			chunk2 <- getNextChunk( reader2, n=100000)
			nNow <- size( chunk1)
			if ( nNow < 100000) hasMore <- FALSE
			if ( nNow < 1) break
	
			readIDs1 <- readID( chunk1)
			readIDs2 <- readID( chunk2)
			# some Sequencing cores will leave mate suffixes...  strip them off
			readIDs1 <- sub( "/1$", "", readIDs1)
			readIDs2 <- sub( "/2$", "", readIDs2)
			refIDs1 <- refID( chunk1)
			refIDs2 <- refID( chunk2)
			seqIDs1 <- refID2seqID( refIDs1, refData=refData)
			seqIDs2 <- refID2seqID( refIDs2, refData=refData)
			seqIDs1[ is.na( seqIDs1)] <- "NoHit"
			seqIDs2[ is.na( seqIDs2)] <- "NoHit"
			pos1 <- position( chunk1)
			pos2 <- position( chunk2)
			isNoHit1 <- which( seqIDs1 == "NoHit")
			if ( length( isNoHit1)) {
				seqDNA1 <- readSeq(chunk1)
				# check for empty adaptor ends
				scores <- pairwiseAlignment( seqDNA1[isNoHit1], rcadapter2, type='local', scoreOnly=TRUE)
				reLabel <- isNoHit1[ scores >= adapterMinScore]
				if (length(reLabel)) seqIDs1[ reLabel] <- "EmptyAdapter"
				if (trackNoHits) {
					smlTbl <- sort( table( seqDNA1[ isNoHit1]), decreasing=T)
					Nkeep <- min( N_NOHIT_KEEP, length(smlTbl))
					nohitTbl1 <- if (is.null(nohitTbl1)) smlTbl[ 1:Nkeep] else mergeTables( nohitTbl1, smlTbl[1:Nkeep])
				}
				# check for Poly N
				isNoHit1 <- setdiff( isNoHit1, reLabel)
				bases <- strsplit( seqDNA1[isNoHit1], split="")
				nPolyN <- sapply( bases, function(x) sum( x == "N"))
				isPolyN <- which( nPolyN > 5)
				if (length(isPolyN)) seqIDs1[ isNoHit1[isPolyN]] <- "PolyN"
			}
			isNoHit2 <- which( seqIDs2 == "NoHit")
			if ( length( isNoHit2)) {
				seqDNA2 <- readSeq(chunk2)
				# check for empty adaptor ends
				scores <- pairwiseAlignment( seqDNA2[isNoHit2], rcadapter1, type='local', scoreOnly=TRUE)
				reLabel <- isNoHit2[ scores >= adapterMinScore]
				if (length(reLabel)) seqIDs2[ reLabel] <- "EmptyAdapter"
				if (trackNoHits) {
					smlTbl <- sort( table( seqDNA2[ isNoHit2]), decreasing=T)
					Nkeep <- min( N_NOHIT_KEEP, length(smlTbl))
					nohitTbl2 <- if (is.null(nohitTbl2)) smlTbl[ 1:Nkeep] else mergeTables( nohitTbl2, smlTbl[1:Nkeep])
				}
				# check for Poly N
				isNoHit2 <- setdiff( isNoHit2, reLabel)
				bases <- strsplit( seqDNA2[isNoHit2], split="")
				nPolyN <- sapply( bases, function(x) sum( x == "N"))
				isPolyN <- which( nPolyN > 5)
				if (length(isPolyN)) seqIDs2[ isNoHit2[isPolyN]] <- "PolyN"
			}
	
			# force the mate pairs to match by their IDs
			if ( ! all( readIDs1 == readIDs2)) {
				both <- intersect( readIDs1, readIDs2)
				wh1 <- match( both, readIDs1)
				wh2 <- match( both, readIDs2)
				refIDs1 <- refIDs1[wh1]
				refIDs2 <- refIDs2[wh2]
				seqIDs1 <- seqIDs1[wh1]
				seqIDs2 <- seqIDs2[wh2]
				pos1 <- pos1[wh1]
				pos2 <- pos2[wh2]
				readIDs1 <- readIDs2 <- both
			}
			# sanity check that we have intact pairs
			if ( length( readIDs1) < (nNow/2)) stop( "Error:  Problem pairing up Mate Read IDs.")
			cat( ".")
			big1 <- c( big1, seqIDs1)
			big2 <- c( big2, seqIDs2)
			bigID <- c( bigID, readIDs1)
			big1pos <- c( big1pos, pos1)
			big2pos <- c( big2pos, pos2)
		}
	
		bamClose( reader1)
		bamClose( reader2)
		cat( "Done.\n")
			
		# each TFKO entry should be two calls of the same Gene
		Npairs <- length(big1)
		cat( "\nAssess pairing success..")
		text1 <- text2 <- rep.int( "GENE", Npairs)
		text1[ grep( "NoHit", big1, fixed=T)] <- "NOHIT"
		text2[ grep( "NoHit", big2, fixed=T)] <- "NOHIT"
		text1[ grep( "EmptyAdapter", big1, fixed=T)] <- "EMPTY_ADAPT"
		text2[ grep( "EmptyAdapter", big2, fixed=T)] <- "EMPTY_ADAPT"
		text1[ grep( "PolyN", big1, fixed=T)] <- "POLY_N"
		text2[ grep( "PolyN", big2, fixed=T)] <- "POLY_N"
		bigStr <- paste( text1, text2, sep=".")
		# now that the various 'No Hit' reasons are set, call those that actually are what we want
		toCall <- which( bigStr == "GENE.GENE")
		isSameGene <- ( big1[toCall] == big2[toCall])
		bigStr[toCall[ isSameGene]] <- "SAME.GENE"
		bigStr[toCall[ ! isSameGene]] <- "DIFFERENT.GENES"

		pairDF <- data.frame( "ReadID"=bigID, "Mate1"=big1, "Mate2"=big2, "Status"=bigStr, 
				"Pos1"=big1pos, "Pos2"=big2pos, stringsAsFactors=FALSE)
		save( pairDF, file=pairFile)
		cat( "\nWrote Read Pair results:  ", sid, "\tN_Pairs: ", Npairs)

		if (trackNoHits && !is.null( nohitTbl1)) {
			nohitFile <- file.path( pairPath, paste( sid, "TopNoHits.Mate1.fasta", sep="."))
			nohitTbl1 <- sort( nohitTbl1, decreasing=T)
			if ( length(nohitTbl1) > N_NOHIT_KEEP) nohitTbl1 <- nohitTbl1[ 1:N_NOHIT_KEEP]
			NFA <- length(nohitTbl1)
			writeFasta( as.Fasta( desc=paste( "Mate1_", 1:NFA, "_Count=", nohitTbl1, sep=""),
					seq=names(nohitTbl1)), nohitFile, line.width=200)
		}
		if (trackNoHits && !is.null( nohitTbl2)) {
			nohitFile <- file.path( pairPath, paste( sid, "TopNoHits.Mate2.fasta", sep="."))
			nohitTbl2 <- sort( nohitTbl2, decreasing=T)
			if ( length(nohitTbl2) > N_NOHIT_KEEP) nohitTbl2 <- nohitTbl2[ 1:N_NOHIT_KEEP]
			NFA <- length(nohitTbl2)
			writeFasta( as.Fasta( desc=paste( "Mate2_", 1:NFA, "_Count=", nohitTbl2, sep=""),
					seq=names(nohitTbl2)), nohitFile, line.width=200)
		}

	} else {
		cat( "\nSkip..  Read Pair results exist already: ", sid)
	}
}


# function to summarize the TFOE pairings, and make a pie chart
summarizeMutantAlignments.TFOE <- function( sid, rowptr=0, results.path=".", makePie=TRUE,
					bowtieMode="Bowtie2") {

	# where the pair results will be found
	pairPath <- file.path( results.path, paste( "PAIRresults", bowtieMode, sep="."))
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))
	if ( ! file.exists( pairFile)) {
		cat( "\niRead Pair File not found: ", sid, "\t", pairFile)
		return( NULL)
	}
	cat( "\nSummarizing Read Pairs for:  ", sid)

	# read those pair results in
	load( pairFile)
	Npairs <- nrow(pairDF)

	# get all the various groups and their counts
	statusTypes <- c( "ANCHOR", "EMPTY_ADAPT", "GENE", "NOHIT", "POLY_N")
	Ntypes <- length( statusTypes)
	statusLevels <- paste( rep( statusTypes, each=Ntypes), rep( statusTypes, times=Ntypes), sep=".")
	statusFac <- factor( pairDF$Status, levels=statusLevels)
	statusSets <- tapply( 1:Npairs, statusFac, function(x) return(x))
	isDoubleAnch <- statusSets[[1]]
	isAnchEmpty <- statusSets[[2]]
	isAnchGene <- statusSets[[3]]
	isAnchNo <- statusSets[[4]]
	isAnchPolyn <- statusSets[[5]]
	isEmptyAnch <- statusSets[[6]]
	isDoubleEmpty <- statusSets[[7]]
	isEmptyGene <- statusSets[[8]]
	isEmptyNo <- statusSets[[9]]
	isEmptyPolyn <- statusSets[[10]]
	isGeneAnch <- statusSets[[11]]
	isGeneEmpty <- statusSets[[12]]
	isDoubleGene <- statusSets[[13]]
	isGeneNo <- statusSets[[14]]
	isGenePolyn <- statusSets[[15]]
	isNoAnch <- statusSets[[16]]
	isNoEmpty <- statusSets[[17]]
	isNoGene <- statusSets[[18]]
	isDoubleNo <- statusSets[[19]]
	isNoPolyn <- statusSets[[20]]
	isPolynAnch <- statusSets[[21]]
	isPolynEmpty <- statusSets[[22]]
	isPolynGene <- statusSets[[23]]
	isPolynNo <- statusSets[[24]]
	isDoublePolyn <- statusSets[[25]]

	nAnch1 <- length( isAnch1 <- c( isDoubleAnch, isAnchEmpty, isAnchGene, isAnchNo, isAnchPolyn))
	nGene1 <- length( isGene1 <- c( isGeneAnch, isGeneEmpty, isDoubleGene, isGeneNo, isGenePolyn))
	nEmpty1 <- length( isEmpty1 <- c( isEmptyAnch, isDoubleEmpty, isEmptyGene, isEmptyNo, isEmptyPolyn))
	nNoHit1 <- length( isNoHit1 <- c( isNoAnch, isNoEmpty, isNoGene, isDoubleNo, isNoPolyn))
	nPolyn1 <- length( isPolyn1 <- c( isPolynAnch, isPolynEmpty, isPolynGene, isPolynNo, isDoublePolyn))
	nAnch2 <- length( isAnch2 <- c( isDoubleAnch, isEmptyAnch, isGeneAnch, isNoAnch, isPolynAnch))
	nGene2 <- length( isGene2 <- c( isAnchGene, isEmptyGene, isDoubleGene, isNoGene, isPolynGene))
	nEmpty2 <- length( isEmpty2 <- c( isAnchEmpty, isDoubleEmpty, isGeneEmpty, isNoEmpty, isPolynEmpty))
	nNoHit2 <- length( isNoHit2 <- c( isAnchNo, isEmptyNo, isGeneNo, isDoubleNo, isPolynNo))
	nPolyn2 <- length( isPolyn2 <- c( isAnchPolyn, isEmptyPolyn, isGenePolyn, isNoPolyn, isDoublePolyn))

	# make a visual of the distribution
	if (makePie) {
		cnts <- c( length(isGeneAnch), length(isAnchGene), length(isAnchNo), length( isNoAnch), length(isGeneNo), length( isNoGene),
				length(isDoubleAnch), length(isDoubleGene), length(isDoubleNo), length( isDoubleEmpty), 
				length(isAnchEmpty), length( isEmptyAnch), length( isGeneEmpty), length( isEmptyGene),
				length( isEmptyNo), length( isNoEmpty), length(isAnchPolyn), length(isPolynAnch), length( isGenePolyn),
				length( isPolynGene), length(isEmptyPolyn), length( isPolynEmpty), length( isNoPolyn),
				length( isPolynNo), length( isDoublePolyn))
		nams <- c( "Gene - Anchor", "Anchor - Gene", "Anchor - NoHit", "NoHit - Anchor", "Gene - NoHit", "NoHit - Gene", 
				"Anchor - Anchor", "Gene - Gene", "NoHit - NoHit", "EmptyAdapt - EmptyAdapt",
				"Anchor - EmptyAdapt", "EmptyAdapt - Anchor", "Gene - EmptyAdapt", "EmptyAdapt - Gene",
				"EmptyAdapt - NoHit", "NoHit - EmptyAdapt", "Anchor - PolyN", "PolyN - Anchor", "Gene - PolyN",
				"PolyN - Gene", "EmptyAdapt - PolyN", "PolyN - EmptyAdapt", "NoHit - PolyN", "PolyN - NoHit",
				"PolyN - PolyN")
		pcts <- cnts * 100 / sum(cnts)
		strPcts <- as.percent( cnts, big.value=sum(cnts))
		names( cnts) <- paste( nams, "  (", strPcts, ")", sep="")
		cols <- c( 3, 3, "gold", "gold", "orange", "orange", 
				"yellow", "pink", 2, 'purple', 
				'blue', "blue", 'dodgerblue', 'dodgerblue', 'slateblue', 'slateblue',
				rep.int( 'sandybrown',4), rep.int( 'saddlebrown',4), rep.int( 'brown',3))
		toShow <- which( pcts >= 0.75)

		pie( cnts[toShow], col=cols[toShow], clock=F, 
				main=paste( "Read Pair Profile:   ", sid, "\nTotal Reads:   ", 
				formatC( Npairs, format="d", big.mark=",")))

		# try to show each mates overall status as color bars along lower edge
		usrLims <- par( 'usr')
		xlo <- usrLims[1]
		xhi <- usrLims[2]
		xwid <- xhi - xlo
		ylo <- usrLims[3]
		yhi <- usrLims[4]
		ywid <- yhi - ylo
		ytop <- ylo + ywid * 0.05
		ymid <- (ylo + ytop) / 2
		xright <- xlo + xwid * 0.33
		rect( xlo, ylo, xright, ytop, border=1, col='red')
		text( xright, ymid, "NoHit", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1+nGene1+nEmpty1+nPolyn1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='brown')
		text( xright, ymid, "PolyN", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1+nGene1+nEmpty1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='purple')
		text( xright, ymid, "Adaptor", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1+nGene1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='green')
		text( xright, ymid, "Gene", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nAnch1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='yellowgreen')
		text( xright, ymid, "Anchor", pos=2, offset=0.2, cex=0.9)
		xleft <- xhi - xwid * 0.33
		rect( xleft, ylo, xhi, ytop, border=1, col='red')
		text( xleft, ymid, "NoHit", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2+nGene2+nEmpty2+nPolyn2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='brown')
		text( xleft, ymid, "PolyN", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2+nGene2+nEmpty2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='purple')
		text( xleft, ymid, "Adaptor", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2+nGene2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='green')
		text( xleft, ymid, "Gene", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nAnch2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='yellowgreen')
		text( xleft, ymid, "Anchor", pos=4, offset=0.2, cex=0.9)
		text( xlo+xwid*0.16, ytop, "Mate 1 Status", pos=3, cex=1, font=2)
		text( xhi-xwid*0.16, ytop, "Mate 2 Status", pos=3, cex=1, font=2)
		subText <- paste( "Valid Read Pairs:   ", formatC( length( c( isGeneAnch, isAnchGene)), 
				format="d", big.mark=","))
		text( 0, ylo, subText, font=2, cex=1.2, pos=3)

		piePath <- file.path( results.path, paste( "PIEplots", bowtieMode, sep="."))
		if ( ! file.exists( piePath)) dir.create( piePath, recursive=T)
		pieFile <- file.path( piePath, paste( sid, "MatePair.Profile.png", sep="."))
		dev.print( png, pieFile, width=900, height=600)
	}

	# final result is the counts of the GeneIDs for the good rows
	geneNames <- c( pairDF$Mate1[ isGeneAnch], pairDF$Mate2[ isAnchGene])
	out <- table( geneNames)
	cat( "\nN_BAM mate pairs considered:             ", Npairs)
	cat( "\nTotal Good 'Anchor - Gene' mate pairs:   ", sum(out))
	cat( "\nTotal Good mate pairs as percentage:    ", as.percent( sum(out), big.value=Npairs))
	cat( "\nN_Detected Genes:                        ", length(out))

	# save some overall metrics
	valid1 <- nAnch1 + nGene1
	valid2 <- nAnch2 + nGene2
	if ( rowptr > 0 && exists( "outStats")) {
		outStats[ rowptr, 'StartingReads'] <<- Npairs
		outStats[ rowptr, 'ValidMate1'] <<- valid1
		outStats[ rowptr, 'Pct_ValidMate1'] <<- round( valid1 * 100 / Npairs, digits=2)
		outStats[ rowptr, 'ValidMate2'] <<- valid2
		outStats[ rowptr, 'Pct_ValidMate2'] <<- round( valid2 * 100 / Npairs, digits=2)
		outStats[ rowptr, 'ValidPairs'] <<- sum(out)
		outStats[ rowptr, 'Pct_ValidPairs'] <<- round( sum(out) * 100 / Npairs, digits=2)
		outStats[ rowptr, 'N_Genes'] <<- length(out)
	}
	return( out)
}


# function to summarize the TFKO pairings, and make a pie chart
summarizeMutantAlignments.TFKO <- function( sid, rowptr=0, results.path=".", makePie=TRUE,
					bowtieMode="Bowtie2") {

	# where the pair results will be found
	pairPath <- file.path( results.path, paste( "PAIRresults", bowtieMode, sep="."))
	pairFile <- file.path( pairPath, paste( sid, "ReadPairs.rda", sep="."))
	if ( ! file.exists( pairFile)) {
		cat( "\nRead Pair File not found: ", sid, "\t", pairFile)
		return( NULL)
	}
	cat( "\nSummarizing Read Pairs for:  ", sid)

	# read those pair results in
	load( pairFile)
	Npairs <- nrow(pairDF)
	if ( ! Npairs) {
		cat( "\nRead Pair File is empty: ", sid, "\t", pairFile)
		return( NULL)
	}

	# get all the various groups and their counts
	statusTypes <- c( "EMPTY_ADAPT", "GENE", "NOHIT", "POLY_N")
	Ntypes <- length( statusTypes)
	statusLevels <- c( "SAME.GENE", "DIFFERENT.GENES", paste( rep( statusTypes, each=Ntypes), rep( statusTypes, times=Ntypes), sep="."))
	#cat( "\nDebug: ", length(statusLevels), statusLevels, "\n")
	statusFac <- factor( pairDF$Status, levels=statusLevels)
	statusSets <- tapply( 1:Npairs, statusFac, function(x) return(x))
	isSameGene <- statusSets[[1]]
	isDiffGene <- statusSets[[2]]
	isDoubleEmpty <- statusSets[[3]]
	isEmptyGene <- statusSets[[4]]
	isEmptyNo <- statusSets[[5]]
	isEmptyPolyn <- statusSets[[6]]
	isGeneEmpty <- statusSets[[7]]
	isGeneGene <- statusSets[[8]]  # should be empty set
	isGeneNo <- statusSets[[9]]
	isGenePolyn <- statusSets[[10]]
	isNoEmpty <- statusSets[[11]]
	isNoGene <- statusSets[[12]]
	isDoubleNo <- statusSets[[13]]
	isNoPolyn <- statusSets[[14]]
	isPolynEmpty <- statusSets[[15]]
	isPolynGene <- statusSets[[16]]
	isPolynNo <- statusSets[[17]]
	isDoublePolyn <- statusSets[[18]]

	nGene1 <- length( isGene1 <- c( isGeneEmpty, isDiffGene, isSameGene, isGeneNo, isGenePolyn))
	nEmpty1 <- length( isEmpty1 <- c( isDoubleEmpty, isEmptyGene, isEmptyNo, isEmptyPolyn))
	nNoHit1 <- length( isNoHit1 <- c( isNoEmpty, isNoGene, isDoubleNo, isNoPolyn))
	nPolyn1 <- length( isPolyn1 <- c( isPolynEmpty, isPolynGene, isPolynNo, isDoublePolyn))
	nGene2 <- length( isGene2 <- c( isEmptyGene, isDiffGene, isSameGene, isNoGene, isPolynGene))
	nEmpty2 <- length( isEmpty2 <- c( isDoubleEmpty, isGeneEmpty, isNoEmpty, isPolynEmpty))
	nNoHit2 <- length( isNoHit2 <- c( isEmptyNo, isGeneNo, isDoubleNo, isPolynNo))
	nPolyn2 <- length( isPolyn2 <- c( isEmptyPolyn, isGenePolyn, isNoPolyn, isDoublePolyn))

	# make a visual of the distribution
	if (makePie) {
		cnts <- c( length( isDiffGene), length( isSameGene), length(isGeneNo), length( isNoGene),
				length(isDoubleNo), length( isDoubleEmpty), length( isGeneEmpty), length( isEmptyGene),
				length( isEmptyNo), length( isNoEmpty), length( isGenePolyn), length( isPolynGene), 
				length(isEmptyPolyn), length( isPolynEmpty), length( isNoPolyn), length( isPolynNo), 
				length( isDoublePolyn))
		nams <- c( "Different Genes", "Same Gene", "Gene - NoHit", "NoHit - Gene", 
				"NoHit - NoHit", "EmptyAdapt - EmptyAdapt", "Gene - EmptyAdapt", "EmptyAdapt - Gene",
				"EmptyAdapt - NoHit", "NoHit - EmptyAdapt", "Gene - PolyN", "PolyN - Gene", 
				"EmptyAdapt - PolyN", "PolyN - EmptyAdapt", "NoHit - PolyN", "PolyN - NoHit",
				"PolyN - PolyN")
		pcts <- cnts * 100 / sum(cnts)
		strPcts <- as.percent( cnts, big.value=sum(cnts))
		names( cnts) <- paste( nams, "  (", strPcts, ")", sep="")
		cols <- c( "pink", 3, "orange", "orange", 
				2, 'purple', "violet", "Violet",
				'slateblue', 'slateblue',
				rep.int( 'saddlebrown',4), rep.int( 'brown',3))
		toShow <- which( pcts >= 0.75)

		pie( cnts[toShow], col=cols[toShow], clock=F, 
				main=paste( "Read Pair Profile:   ", sid, "\nTotal Reads:   ", 
				formatC( Npairs, format="d", big.mark=",")))

		# try to show each mates overall status as color bars along lower edge
		usrLims <- par( 'usr')
		xlo <- usrLims[1]
		xhi <- usrLims[2]
		xwid <- xhi - xlo
		ylo <- usrLims[3]
		yhi <- usrLims[4]
		ywid <- yhi - ylo
		ytop <- ylo + ywid * 0.05
		ymid <- (ylo + ytop) / 2
		xright <- xlo + xwid * 0.33
		rect( xlo, ylo, xright, ytop, border=1, col='red')
		text( xright, ymid, "NoHit", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nGene1+nEmpty1+nPolyn1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='brown')
		text( xright, ymid, "PolyN", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nGene1+nEmpty1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='purple')
		text( xright, ymid, "Adaptor", pos=2, offset=0.2, cex=0.9)
		xright <- xlo + (xwid*0.33) * ((nGene1)/Npairs)
		rect( xlo, ylo, xright, ytop, border=1, col='green')
		text( xright, ymid, "Gene", pos=2, offset=0.2, cex=0.9)
		xleft <- xhi - xwid * 0.33
		rect( xleft, ylo, xhi, ytop, border=1, col='red')
		text( xleft, ymid, "NoHit", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nGene2+nEmpty2+nPolyn2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='brown')
		text( xleft, ymid, "PolyN", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nGene2+nEmpty2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='purple')
		text( xleft, ymid, "Adaptor", pos=4, offset=0.2, cex=0.9)
		xleft <- xhi - (xwid*0.33) * ((nGene2)/Npairs)
		rect( xleft, ylo, xhi, ytop, border=1, col='green')
		text( xleft, ymid, "Gene", pos=4, offset=0.2, cex=0.9)
		text( xlo+xwid*0.16, ytop, "Mate 1 Status", pos=3, cex=1, font=2)
		text( xhi-xwid*0.16, ytop, "Mate 2 Status", pos=3, cex=1, font=2)
		subText <- paste( "Valid Read Pairs:   ", formatC( length( isSameGene), 
				format="d", big.mark=","))
		text( 0, ylo, subText, font=2, cex=1.2, pos=3)

		piePath <- file.path( results.path, paste( "PIEplots", bowtieMode, sep="."))
		if ( ! file.exists( piePath)) dir.create( piePath, recursive=T)
		pieFile <- file.path( piePath, paste( sid, "MatePair.Profile.png", sep="."))
		dev.print( png, pieFile, width=900, height=600)
	}

	# final result is the counts of the 'SameGene' for the good rows
	geneNames <- pairDF$Mate1[ isSameGene]
	out <- table( geneNames)
	cat( "\nN_BAM mate pairs considered:             ", Npairs)
	cat( "\nTotal Good 'Same Gene' mate pairs:   ", sum(out))
	cat( "\nTotal Good mate pairs as percentage:    ", as.percent( sum(out), big.value=Npairs))
	cat( "\nN_Detected Genes:                        ", length(out))

	# save some overall metrics
	valid1 <- nGene1
	valid2 <- nGene2
	if ( rowptr > 0 && exists( "outStats")) {
		outStats[ rowptr, 'StartingReads'] <<- Npairs
		outStats[ rowptr, 'ValidMate1'] <<- valid1
		outStats[ rowptr, 'Pct_ValidMate1'] <<- round( valid1 * 100 / Npairs, digits=2)
		outStats[ rowptr, 'ValidMate2'] <<- valid2
		outStats[ rowptr, 'Pct_ValidMate2'] <<- round( valid2 * 100 / Npairs, digits=2)
		outStats[ rowptr, 'ValidPairs'] <<- sum(out)
		outStats[ rowptr, 'Pct_ValidPairs'] <<- round( sum(out) * 100 / Npairs, digits=2)
		outStats[ rowptr, 'N_Genes'] <<- length(out)
	}
	return( out)
}


writeResultTables <- function( tbl, samples, results.path, nExpectedGenes, poolGenes, experiment, 
				bowtieMode, makeROC=FALSE, runName="") {

	# turn these raw counts of all genes into log2 RPMHEG, etc.
	outRPMlog2 <- tbl
	expectRPMlog2 <- expect <- tbl[ which( rownames(tbl) %in% poolGenes), ]

	cat( "\nNormalizing to RPM_HEG units..\n")
	nSamples <- ncol(tbl)
	nGenes <- nrow(tbl)
	if ( nrow(samples) != nSamples) stop( "Sample / matrix sizing error!")

	for ( i in 1:nSamples) {
		v <- tbl[ , i]
		expectedGenesFactor <- nExpectedGenes[i] / 100
		vnorm <- v * expectedGenesFactor * 1000000 / sum(v)
		outRPMlog2[ , i] <- log2( vnorm + 1)
		v <- expect[ , i]
		vnorm <- v * expectedGenesFactor * 1000000 / sum(v)
		expectRPMlog2[ , i] <- log2( vnorm + 1)
		cat( "\r", colnames(tbl)[i])
	}

	outfile <- file.path( results.path, paste( "AllGenes.RawCounts", experiment, bowtieMode, "txt", sep="."))
	out <- data.frame( "GeneID"=rownames(tbl), tbl, stringsAsFactors=F)
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.RawCounts", experiment, bowtieMode, "txt", sep="."))
	expect <- data.frame( "GeneID"=rownames(expect), expect, stringsAsFactors=F)
	write.table( expect, expectfile, sep="\t", quote=F, row.names=F)

	if ( makeROC) {
		makeROCimage( outRPMlog2, poolGenes=poolGenes, text="Log2 RPMHEG", experiment=experiment, 
				bowtieMode=bowtieMode, results.path=results.path)
	}
	outfile <- file.path( results.path, paste( "AllGenes.Log2RPMHEG", experiment, bowtieMode, "txt", sep="."))
	out <- data.frame( "GeneID"=rownames(outRPMlog2), outRPMlog2, stringsAsFactors=F)
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	expectfile <- file.path( results.path, paste( "PoolGenes.Log2RPMHEG", experiment, bowtieMode, "txt", sep="."))
	expect <- data.frame( "GeneID"=rownames(expectRPMlog2), expectRPMlog2, stringsAsFactors=F)
	write.table( expect, expectfile, sep="\t", quote=F, row.names=F)

	# lastly, make a table with all meta data and all genes for modeling etc.
	completeGeneTbl <- read.delim( CompleteGeneFile, as.is=T)
	completeGenes <- sort( unique( completeGeneTbl$GeneID))
	outCounts <- matrix( NA, nrow=ncol(expectRPMlog2), ncol=length(completeGenes))
	colnames(outCounts) <- completeGenes
	rownames(outCounts) <- samples$SampleID
	# first the log2 RPMHEG counts
	where <- match( rownames(expectRPMlog2), completeGenes)
	if ( any( is.na(where))) {
		cat( "\n\nSome GeneIDs not in the 'CompleteGene' set:\n")
		cat( "\nN: ", sum(is.na(where)), "  Who:  ", rownames(expectRPMlog2)[is.na(where)])
	}
	for (i in 1:nGenes) {
		if ( is.na( where[i])) next
		outCounts[ , where[i]] <- expectRPMlog2[ i, ]
	}
	# then all the meta data (SampleID through Mutant Pool)
	outMeta <- samples[ , MetaDataColumns]

	# Tige wants the run name and the raw read counts
	readCnts <- apply( tbl, MARGIN=2, sum, na.rm=T)
	out <- data.frame( "SampleID"=outMeta[ ,1], "Run"=runName, outMeta[ ,2:ncol(outMeta)], "ReadCount"=readCnts, 
				round( outCounts, digits=4), stringsAsFactors=FALSE)
	rownames(out) <- 1:nrow(out)
	outfile <- file.path( results.path, paste( "FinalResults", experiment, "txt", sep="."))
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
}


combineResults <- function( folder=SampleKeyFolder, mode=c("TFOE","TFKO")) {

	# get all the sample key files
	keyFiles <- dir( folder, patt="^SamplesKey.+\\.txt$", full.name=T)
	nFiles <- length(keyFiles)
	cat( "\nN_Sample Key files: ", nFiles)

	out <- data.frame()
	NC.out <- 0
	for ( i in 1:nFiles) {
		samples <- read.delim( keyFiles[i], as.is=T)
		runName <- extractRunName( keyFiles[i])
		results.path <- createResultsFolder( keyFiles[i], create=FALSE)
		if ( is.null( results.path)) {
			cat( "\nResults folder not found: ", runName)
			next
		}
		experiments <- sort( unique( samples$Experiment))
		for (e in experiments) {
			f <- file.path( results.path, paste( "FinalResults", e, "txt", sep="."))
			if ( ! file.exists(f)) {
				cat( "\nFinal Results File not found: ", f)
				next
			}
			thisDF <- read.delim( f, as.is=T)

			# check before blindly combining
			if ( NC.out == 0) {
				NC.out <- ncol(thisDF)
			} else {
				if ( ncol(thisDF) != NC.out) {
					cat( "\nWrong size results file..  Skipping: ", basename(f))
					next
				}
			}
			out <- rbind( out, thisDF)
			cat( "\nRun: ", runName, "  Experiment: ", e, "  N_Samples: ", nrow(thisDF))
		}
	}
	nExp <- length( unique( out$Experiment))
	nSamples <- nrow(out)
	cat( "\nN_Experiments: ", nExp)
	cat( "\nN_Samples:     ", nSamples)

	outfile <- file.path( ResultsFolder, "FinalResults.All.TFOE.Log2RPMHEG.txt")
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
}


getSampleLabel <- function( samples, mode=c("isolate", "group", "ID")) {

	mode <- match.arg( mode)

	# using the named columns, not the SampleID anymore
	N <- nrow(samples)
	pool <- samples$MutantPool
	nPools <- length( unique( pool))
	class <- samples$Class
	nClass <- length( unique( class))
	sub1 <- samples$SubClass1
	nSub1 <- length( unique( sub1[sub1 != ""]))
	sub2 <- samples$SubClass2
	# sub clase 2 may have decimal quantities
	sub2 <- sub( "0p", "0.", sub2, fixed=T)
	nSub2 <- length( unique( sub2[sub2 != ""]))

	# try to build a meaningful key, but not full of empty terms
	key <- rep.int( "", N)
	if (nPools > 1) key <- pool
	seps <- ifelse( nchar(key) > 0, "_", "")
	if (nClass > 1) key <- paste( key, seps, class, sep="")
	# if we have no details yet, use the class
	if ( all( key == "")) key <- class
	seps <- ifelse( nchar(key) > 0, "_", "")
	if (nSub1 > 0) key <- paste( key, seps, sub1, sep="")
	seps <- ifelse( nchar(key) > 0, "_", "")
	if (nSub2 > 0) key <- paste( key, seps, sub2, sep="")
	sampleKey <- key
	if ( mode == "group") return( sampleKey)

	induced <- ifelse( samples$Induced, "Ind", "Un")
	day <- paste( "D", samples$Day, sep="")
	repl <- paste( "r", samples$Replicate, sep="")
	nRepl <- length( unique( repl))
	txt <- paste( day, induced, sep="_")
	if ( nRepl > 1) txt <- paste( txt, repl, sep="_")
	# prepend sub class info to the text
	seps <- ifelse( nchar(sub2) > 0, "_", "")
	if (nSub2 > 1) txt <- paste( sub2, seps, txt, sep="")
	seps <- ifelse( nchar(sub1) > 0, "_", "")
	if (nSub1 > 1) txt <- paste( sub1, seps, txt, sep="")
	sampleText <- txt
	if ( mode == "isolate") return( sampleText)

	return( samples$SampleID)
}


plotAlignmentOverview <- function( m=outStats, samples, bowtieMode="Bowtie2", results.path=".") {

	require( plotrix)
	par( mai=c( 1, 1, 0.8, 0.4))

	x <- m[ ,'N_Genes']
	y <- m[,'Pct_ValidPairs']
	nPts <- nrow(m)
	sampleKey <- getSampleLabel( samples, mode="group")
	sampleText <- getSampleLabel( samples, mode="isolate")
	keys <- sort( unique( sampleKey))
	nKeys <- length( keys)
	colorSet <- rainbow( nKeys, end=0.8)
	col <- colorSet[ match( sampleKey, keys)]
	
	# set up, and size by how many to draw
	xlim <- range( x) * c( 0.3, 1.25)
	txt.cex <- 1
	if ( nPts > 6) txt.cex <- 0.9
	if ( nPts > 12) txt.cex <- 0.87
	if ( nPts > 18) txt.cex <- 0.84
	if ( nPts > 24) txt.cex <- 0.81
	if ( nPts > 30) txt.cex <- 0.78
	if ( nPts > 40) txt.cex <- 0.75
	if ( nPts > 60) txt.cex <- 0.72
	leg.cex <- 1
	if ( nKeys > 4) leg.cex <- 0.9
	if ( nKeys > 8) leg.cex <- 0.82
	if ( nKeys > 12) leg.cex <- 0.74
	if ( nKeys > 16) leg.cex <- 0.66

	plot( x, y, xlab="Number of Genes Detected per Sample", xlim=xlim, 
		ylab="Pct Valid Mate Pairs per Sample", ylim=c(0,100),
		main=paste( "Alignment Success Overview:    ", basename(results.path),
		"\n(program used:   '", bowtieMode, "')", sep=""),
		pch=21, cex=1.5, bg=col, col=1, font.lab=2, font.axis=2)

	# replot in random order to show colors better
	ord <- sample( length(x))
	points( x[ord], y[ord], pch=21, cex=1.5, bg=col[ord])
	thigmophobe.labels( x, y, sampleText, cex=txt.cex, font=1)

	legend( 'bottomright', keys, pch=21, pt.bg=colorSet, cex=leg.cex, pt.cex=1.3)
	outfile <- file.path( results.path, paste( "AlignmentSuccessOverview", bowtieMode, "png", sep="."))
	dev.print( png, outfile, width=900, height=600)
}


replotAlignmentOverview <- function( sampleKeyFile="SamplesKey.Run1.txt", bowtieMode=c( "Bowtie2","Bowtie1")) {

	# pull back in the results summary and plot it
	samples <- loadSampleKeyFile( sampleKeyFile)
	results.path <- createResultsFolder( sampleKeyFile)
	bowtieMode <- match.arg( bowtieMode)

	# get the file of stats
	statfile <- file.path( results.path, paste( "AlignmentStatistics", bowtieMode, "txt", sep="."))
	outStats <- read.delim( statfile, as.is=T)

	# call the plotter
	plotAlignmentOverview( outStats, samples, bowtieMode, results.path=results.path)
}


makeROCimage <- function( out, poolGenes, text, experiment, bowtieMode, results.path=".") {

	cat( "\nMaking ROC plot..")
	isPool <- which( rownames(out) %in% poolGenes)
	isBad <- setdiff( 1:nrow(out), isPool)

	goodScores <- as.vector( out[ isPool, ])
	badScores <- as.vector( out[ isBad, ])

	textLabel <- paste( "Experiment:  ", experiment, "      Abundance Units:  ", text)
	ans <- duffy.ROC( goodScores, badScores, label=textLabel)

	plotFile <- file.path( results.path, paste( "ROC", experiment, gsub(" ","",text), "png", sep="."))
	dev.print( png, plotFile, width=800,  height=700)
	return( ans)
}


loadSampleKeyFile <- function( sampleKeyFile) {

	sampleKeyFile <- file.path( SampleKeyFolder, sampleKeyFile)
	if ( ! file.exists( sampleKeyFile)) stop( "'Sample key file' not found..  Tried: ", sampleKeyFile)

	allSamples <- read.delim( sampleKeyFile, as.is=T)

	# force any NA to be blanks
	allSamples <- cleanMetaDataFields( allSamples)
	if ( "FailSeq" %in% colnames(allSamples)) {
		toDrop <- which( toupper(as.character(allSamples$FailSeq)) %in% c( "T", "TRUE"))
		if ( length( toDrop)) allSamples <- allSamples[ -toDrop, ]
	}
	nSamples <- nrow( allSamples)
	if ( nSamples < 1) stop( "No non-excluded samples!")

	if ( ! all( SampleKeyColumns %in% colnames(allSamples))) {
		cat( "\nSample Key file must have the following columns:  ", SampleKeyColumns)
		stop()
	}

	cat( "\n\nSample Key: ", sampleKeyFile, "\nN_Samples:      ", nSamples, "\n")
	cat( "\nExperiment Names: ", sort( unique( allSamples$Experiment)))
	cat( "\nMutant Pools:     ", sort( unique( allSamples$MutantPool)))

	# test the existance of the files
	f1 <- file.path( allSamples$FastqPath, allSamples$File1)
	notfound <- which( ! file.exists( f1))
	if ( length(notfound)) {
		cat( "\nMate 1 files not found: ", length(notfound), "\n", allSamples$File1[notfound])
	}
	f2 <- file.path( allSamples$FastqPath, allSamples$File2)
	notfound <- which( ! file.exists( f2))
	if ( length(notfound)) {
		cat( "\nMate 2 files not found: ", length(notfound), "\n", allSamples$File2[notfound])
	}
	return( allSamples)
}


extractRunName <- function( sampleKeyFile) {

	runName <- sub( "^SamplesKey.", "", basename(sampleKeyFile))
	runName <- sub( ".txt$", "", runName)
	return( runName)
}


createResultsFolder <- function( sampleKeyFile, create=TRUE) {

	runName <- extractRunName( sampleKeyFile)
	results.path <- file.path( "Results", runName)
	if ( ! file.exists( results.path)) {
		if (create) {
			dir.create( results.path, recursive=T)
		} else {
			return(NULL)
		}
	}
	return( results.path)
}


getMutantPoolGenes <- function( myPool) {

	poolFile <- paste( myPool, "Primers.txt", sep=".")
	poolFile <- file.path( "/cidr/primary/sherman/bigdata/HTS.MutantPoolDeconvolution/BobMorrison_MutantPipeline/CustomMutantGenomes", poolFile)
	if ( ! file.exists( poolFile)) stop( paste( "Failed to find Mutant Pool Primers file: ", poolFile))
	poolTable <- read.delim( poolFile, as.is=T)
	if ( ! ("GeneID" %in% colnames(poolTable))) stop( "Mutant Pool Primer file must have a 'GeneID' column")
	poolGenes <- poolTable$GeneID
	return( poolGenes)
}


cleanMetaDataFields <- function( tbl) {

	# force any NA to be blanks
	if ( "FailSeq" %in% colnames(tbl)) tbl$FailSeq[ is.na( tbl$FailSeq)] <- ""
	if ( "Class" %in% colnames(tbl)) tbl$Class[ is.na( tbl$Class)] <- ""
	if ( "SubClass1" %in% colnames(tbl)) tbl$SubClass1[ is.na( tbl$SubClass1)] <- ""
	if ( "SubClass2" %in% colnames(tbl)) tbl$SubClass2[ is.na( tbl$SubClass2)] <- ""
	return( tbl)
}


createClassDescriptor <- function( tbl) {

	# try to use the Class/SubClass data to make a useful descriptor for plots, etc.
	tbl <- cleanMetaDataFields( tbl)
	
	# if one field, use it...
	classes <- setdiff( unique( tbl$Class), "")
	subclass1 <- setdiff( unique( tbl$SubClass1), "")
	subclass2 <- setdiff( unique( tbl$SubClass2), "")
	descriptor <- if (length(classes) == 1) classes[1] else "MultiClass"
	if ( length(subclass1)) {
		suffix <- if ( length(subclass1) == 1) subclass1[1] else "MultiSubClass1"
		descriptor <- paste( descriptor, suffix, sep="_")
	}
	if ( length(subclass2)) {
		suffix <- if ( length(subclass2) == 1) subclass2[1] else "MultiSubClass2"
		descriptor <- paste( descriptor, suffix, sep="_")
	}
	return( descriptor)
}


model.TFOE.GrowthDefects <- function( tbl, min.value=MIN_LOG2_RPMHEG, makePlots=FALSE, 
					plot.path=".", plot.prefix=NULL, asPNG=makePlots) {

	# given a data frame of final results (large table with columns of meta data & Mutant gene log2RPMHEG data)
	# calculate all growth defects from 'Induced' vs 'Uninduced'
	# this routine assumes all rows go into the model!
	# if you want to subset by Class / SubClass / etc., do that in a parent function that calls this.

	if ( ! all( MetaDataColumns %in% colnames(tbl))) {
		cat( "\nMissing some required Meta Data Columns..  Invalid input..")
		stop()
	}
	saveModelInput <<- tbl

	# we will run the same model against each gene column, after dropping out any meta data columns
	geneNames <- setdiff( colnames(tbl), ModelMetaDataColumns)
	nGenes <- length( geneNames)

	# and we will model each Day (past Day 1) on its own, 
	daysToModel <- setdiff( sort( unique( tbl$Day)), DAY_ZERO)
	nDays <- length( daysToModel)

	# storage for data frames and model results in case we draw
	dfByDayIND <- dfByDayUNIND <- vector( length=nDays, mode="list")
	names(dfByDayIND) <- names(dfByDayUNIND) <- daysToModel
	modelAnsByDayIND <- modelAnsByDayUNIND <- vector( length=nDays, mode="list")
	names(modelAnsByDayIND) <- names(modelAnsByDayUNIND) <- daysToModel

	# set up results storage
	growthRate <- growthPvalue <- matrix( NA, nrow=nDays, ncol=nGenes)
	rownames(growthRate) <- rownames(growthPvalue) <- paste( "Day", daysToModel, sep="")
	colnames(growthRate) <- colnames(growthPvalue) <- geneNames

	# by rule, all the 'DayZero' data will be used by both Induced and Uninduced
	dayZeroRows <- which( tbl$Day %in% DAY_ZERO)
	inducedRows <- which( tbl$Induced)
	uninducedRows <- which( ! tbl$Induced)

	# if we will make plots, set up for what we need
	if (makePlots) {
		if ( ! file.exists(plot.path)) dir.create( plot.path, recursive=T)
		if ( is.null( plot.prefix)) {
			plot.prefix <- createClassDescriptor(tbl)
		}
	}

	# ready to loop over all Days/Genes
	for ( j in 1:nGenes) {
		thisGene <- geneNames[j]

		# zero out storage passed to plotter
		for ( i in 1:nDays) {
			dfByDayIND[[i]] <- NULL
			dfByDayUNIND[[i]] <- NULL
			modelAnsByDayIND[[i]] <- NULL
			modelAnsByDayUNIND[[i]] <- NULL
		}
		haveDataToPlot <- FALSE

		for ( i in 1:nDays) {
			thisDay <- daysToModel[i]
			myRows <- which( tbl$Day == thisDay)

			# make the small data frames for both models
			wantColumns <- c( thisGene, "Day", "Class", "SubClass1", "SubClass2")
			wantRowsIND <- c( dayZeroRows, intersect( inducedRows, myRows))
			wantRowsUNIND <- c( dayZeroRows, intersect( uninducedRows, myRows))
			smlIND <- tbl[ wantRowsIND, wantColumns]
			smlUNIND <- tbl[ wantRowsUNIND, wantColumns]
			colnames(smlIND) <- colnames(smlUNIND) <- c( "Value", wantColumns[2:length(wantColumns)])

			# values below a threshold are not real -- set to NA
			smlIND$Value[ smlIND$Value < min.value] <- NA
			smlUNIND$Value[ smlUNIND$Value < min.value] <- NA

			# there may be no data values for some genes
			if ( all( is.na( smlIND$Value[ smlIND$Day == thisDay]))) next

			# there may be no uninduced data for this day, 
			# if so propagate the day zero data to be this day's control
			nThisDay <- sum( smlUNIND$Day ==  thisDay)
			if ( nThisDay == 0) {
				extraData <- smlUNIND
				extraData$Day <- thisDay
				smlUNIND <- rbind( smlUNIND, extraData)
			}
			if ( all( is.na( smlUNIND$Value[ smlUNIND$Day == thisDay]))) next

			# if only a single value, the modeling will bust, extend if needed
			smlIND <- forceEnoughModelData( smlIND)
			if ( is.null(smlIND)) next
			smlUNIND <- forceEnoughModelData( smlUNIND)
			if ( is.null(smlUNIND)) next

			# ready to run the model, save the values for drawing
			dfByDayIND[[i]] <- smlIND
			dfByDayUNIND[[i]] <- smlUNIND
			saveIND <<- smlIND
			saveUNIND <<- smlUNIND

			# call the linear model for both 
			modelAnsByDayIND[[i]] <- ansIND <- lm( Value ~ Day, data=smlIND)
			modelAnsByDayUNIND[[i]] <- ansUNIND <- lm( Value ~ Day, data=smlUNIND)

			# see how different they are, first by LM slope
			diffAns <- lmSlopeDifference( ansUNIND, ansIND)
			lmSlope <- diffAns$difference
			lmPval <- diffAns$p.value
			if ( is.na(lmPval) || is.nan(lmPval)) lmPval <- 1

			# also do simple Student's T test on just the day of interest
			unY <- subset( smlUNIND, Day == thisDay)$Value
			indY <- subset( smlIND, Day == thisDay)$Value
			ttestAns <- t.test( unY, indY)
			ttestSlope <- diff( ttestAns$estimate) / thisDay
			ttestPval <- ttestAns$p.value
			if ( is.na(ttestPval) || is.nan(ttestPval)) ttestPval <- 1
			
			# final answer is the average of the 2 methods
			finalGrowth <- growthRate[i,j] <- mean( c( lmSlope, ttestSlope))
			finalPvalue <- growthPvalue[i,j] <- logmean( c( lmPval, ttestPval))
			haveDataToPlot <- TRUE
			cat( "\r", j, thisGene, thisDay, nrow(smlIND), nrow(smlUNIND), lmSlope, 
					ttestSlope, lmPval, ttestPval)
		}

		# we may do a plot...
		if ( haveDataToPlot && makePlots) {
			plotOneGene.TFOE( thisGene, dfByDayIND, dfByDayUNIND, modelAnsByDayIND, modelAnsByDayUNIND, 
					plot.path=plot.path, plot.prefix=plot.prefix, asPNG=asPNG, 
					showGrowth=finalGrowth, showPvalue=finalPvalue)
		}
	}
	return( list( "GrowthDefect"=growthRate, "P_Value"=growthPvalue))
}


model.TFKO.GrowthDefects <- function( tbl, min.value=MIN_LOG2_RPMHEG, makePlots=FALSE, 
					plot.path=".", plot.prefix=NULL, asPNG=makePlots) {

	# given a data frame of final results (large table with columns of meta data & Mutant gene log2RPMHEG data)
	# calculate all growth defects over time
	# this routine assumes all rows go into the model!
	# if you want to subset by Class / SubClass / etc., do that in a parent function that calls this.

	if ( ! all( MetaDataColumns %in% colnames(tbl))) {
		cat( "\nMissing some required Meta Data Columns..  Invalid input..")
		stop()
	}

	saveModelInput <<- tbl

	# we will run the same model against each gene column
	geneNames <- setdiff( colnames(tbl), ModelMetaDataColumns)
	nGenes <- length( geneNames)

	# and we will model each Day (past Day 1) on its own, 
	daysToModel <- setdiff( sort( unique( tbl$Day)), DAY_ZERO)
	nDays <- length( daysToModel)

	# storage for data frames and model results in case we draw
	dfByDay <- vector( length=nDays, mode="list")
	names(dfByDay) <- daysToModel
	modelAnsByDay <- vector( length=nDays, mode="list")
	names(modelAnsByDay) <- daysToModel

	# set up results storage
	growthRate <- growthPvalue <- matrix( NA, nrow=nDays, ncol=nGenes)
	rownames(growthRate) <- rownames(growthPvalue) <- paste( "Day", daysToModel, sep="")
	colnames(growthRate) <- colnames(growthPvalue) <- geneNames

	# by rule, all the 'DayZero' data will be used as the negative control
	dayZeroRows <- which( tbl$Day %in% DAY_ZERO)

	# if we will make plots, set up for what we need
	if (makePlots) {
		if ( ! file.exists(plot.path)) dir.create( plot.path, recursive=T)
		if ( is.null( plot.prefix)) {
			plot.prefix <- createClassDescriptor(tbl)
		}
	}

	# ready to loop over all Days/Genes
	for ( j in 1:nGenes) {
		thisGene <- geneNames[j]

		# zero out storage passed to plotter
		for ( i in 1:nDays) {
			dfByDay[[i]] <- NULL
			modelAnsByDay[[i]] <- NULL
		}
		haveDataToPlot <- FALSE

		for ( i in 1:nDays) {
			thisDay <- daysToModel[i]
			myRows <- which( tbl$Day == thisDay)

			# make the small data frame for the model
			wantColumns <- c( thisGene, "Day", "Class", "SubClass1", "SubClass2")
			wantRows <- c( dayZeroRows, myRows)
			sml <- tbl[ wantRows, wantColumns]
			colnames(sml) <- c( "Value", wantColumns[2:length(wantColumns)])

			# values below a threshold are not real -- set to NA
			sml$Value[ sml$Value < min.value] <- NA

			# there may be no data values for some genes
			if ( all( is.na( sml$Value[ sml$Day == thisDay]))) next

			# if only a single value, the modeling will bust, extend if needed
			sml <- forceEnoughModelData( sml)
			if ( is.null(sml)) next
			saveSmall <<- sml

			# ready to run the model, save the values for drawing
			dfByDay[[i]] <- sml
			names(dfByDay)[i] <- paste( "Day_",thisDay,sep="")

			# call the linear model
			modelAnsByDay[[i]] <- ans <- lm( Value ~ Day, data=sml)
			names(modelAnsByDay)[i] <- paste( "Day_",thisDay,sep="")
			saveLM <<- ans

			# no comparison between 2 models, just the one slope -- is it flat?
			ansSummary <- summary( ans)
			cofs <- coef( ansSummary)
			lmSlope <- cofs[ 2, 1]
			lmPval <- cofs[2, 4]
			if ( is.na(lmPval) || is.nan(lmPval)) lmPval <- 1

			# also do simple Student's T test on just the day of interest
			Y0 <- subset( sml, Day %in% DAY_ZERO)$Value
			Y <- subset( sml, Day == thisDay)$Value
			ttestAns <- t.test( Y0, Y)
			ttestSlope <- diff( ttestAns$estimate) / thisDay
			ttestPval <- ttestAns$p.value
			if ( is.na(ttestPval) || is.nan(ttestPval)) ttestPval <- 1
			
			# final answer is the average of the 2 methods
			finalGrowth <- growthRate[i,j] <- mean( c( lmSlope, ttestSlope))
			finalPvalue <- growthPvalue[i,j] <- logmean( c( lmPval, ttestPval))
			haveDataToPlot <- TRUE
			cat( "\r", j, thisGene, thisDay, nrow(sml), lmSlope, ttestSlope, lmPval, ttestPval)
		}

		# we may do a plot...
		if ( haveDataToPlot && makePlots) {
			plotOneGene.TFKO( thisGene, dfByDay, modelAnsByDay, 
					plot.path=plot.path, plot.prefix=plot.prefix, asPNG=asPNG, 
					showGrowth=finalGrowth, showPvalue=finalPvalue)
		}
	}
	return( list( "GrowthDefect"=growthRate, "P_Value"=growthPvalue))
}


forceEnoughModelData <- function( df) {

	# make sure every day has at least 2 observations
	dayFac <- factor( df$Day)
	nGood <- tapply( !is.na(df$Value), dayFac, sum)
	if ( all( nGood > 1)) return(df)

	df$DrawIt <- TRUE
	out <- df

	for (i in 1:nlevels(dayFac)) {
		if (nGood[i] > 1) next
		thisDay <- as.numeric( levels(dayFac)[i])
		curVal <- df$Value[ df$Day == thisDay]
		goodVal <- curVal[ !is.na(curVal)][1]
		if ( is.na(goodVal)) return(NULL)
		extraRow <- df[ which( df$Value == goodVal & df$Day == thisDay)[1], ]
		extraRow$Value <- jitter(goodVal, factor=1, amount=0.5)
		extraRow$Day <- thisDay
		extraRow$DrawIt <- FALSE
		out <- rbind( out, extraRow)
	}
	return( out)
}


plotOneGene.TFOE <- function( gene, dfsIND, dfsUNIND, modelsIND, modelsUNIND, 
			plot.path=".", plot.prefix="", asPNG=FALSE, 
			showGrowth=NULL, showPvalue=NULL) {

	checkX11( bg="white", width=10, height=7)
	par( mai=c( 1,1,0.8,0.4))
	fileRoot <- paste( gene, plot.prefix, sep="_")

	# scan all the data frames for plot extents
	allDays <- allValues <- allClass <- vector()
	nModels <- length( dfsIND)
	if (nModels < 1) return()
	for (i in 1:nModels) {
		smlIND <- dfsIND[[i]]
		if ( is.null(smlIND)) next
		smlUNIND <- dfsUNIND[[i]]
		if ( is.null(smlUNIND)) next
		# small chance of extra non-drawn data to prevent breaking the model fitting
		if ( "DrawIt" %in% colnames(smlIND)) smlIND <- subset( smlIND, DrawIt)
		if ( "DrawIt" %in% colnames(smlUNIND)) smlUNIND <- subset( smlUNIND, DrawIt)
		allDays <- c( allDays, smlIND$Day)
		allDays <- c( allDays, smlUNIND$Day)
		allValues <- c( allValues, smlIND$Value)
		allValues <- c( allValues, smlUNIND$Value)
	}
	xlim <- range( allDays, na.rm=TRUE) + c( -0.5,0.5)
	ylim <- range( allValues, na.rm=TRUE)
	ylim[2] <- ylim[2] + diff(ylim)*0.35
	if ( ylim[1] > 0) ylim[1] <- ylim[1] * 0.9

	mainText <- paste( "Mutant TFOE Timecourse:    ", plot.prefix, "\n", gene, " -  ", gene2Product(gene))
	plot( 1, 1, type='n', main=mainText, xlim=xlim, ylim=ylim, xlab="Time  (days)", xaxt='n',
			ylab="Abundance   Log2(RPMHEG)", font.axis=2, font.lab=2)
	axis( side=1, at=sort( unique( allDays)), font=2)


	# local function to draw one "Linear Model" object...
	showOneTimeCourse <- function( df, lmAns, thisCol, thisPCH=21, pt.cex=1, showZero=TRUE, showLine=TRUE) {

		if ( is.null(df)) return()

		thisX <- df$Day
		thisY <- df$Value
		xFac <- factor( thisX)
		yMeans <- rep.int( 0, nlevels(xFac))
		saveDF <<- df
		
		# perhaps do not redraw the zero points (where 'Zero' includes day 1)
		toShow <- 1:length(thisX)
		if ( !showZero) toShow <- which( ! (thisX %in% DAY_ZERO))

		points( jitter(thisX,factor=0.2)[toShow], thisY[toShow], pch=thisPCH, bg=thisCol, cex=pt.cex)
		yMeans <- tapply( thisY, xFac, mean)
		mX <- as.numeric( names(yMeans))
		for ( who in which( yMeans > 0)) {
			if ( !showZero && mX[who] %in% DAY_ZERO) next
			lines( mX[who]+c(-0.2,0.2), rep.int(yMeans[who],2), col=thisCol, lwd=2)
		}

		if ( is.null(lmAns)) return()
		if ( showLine && length(yMeans) > 1) {
			abline.segment( reg=lmAns, x1=min(mX), x2=max(mX), col=thisCol, lty=2, lwd=2)
		}
	}


	# OK, all set up, draw those time courses
	colUN <- "dodgerblue"
	colIND <- "red"
	
	# for the uninduced, only show the day zero and best fit for the last model
	for ( i in 1:nModels) {
		if ( is.null( dfsUNIND[[i]])) next
		showZero <- showLine <- (i == nModels)
		showOneTimeCourse( dfsUNIND[[i]], modelsUNIND[[i]], colUN, showZero=showZero, showLine=showLine)
	}

	# for the induced, never show the day zero and always show the best fit 
	for ( i in 1:nModels) {
		if ( is.null( dfsIND[[i]])) next
		showOneTimeCourse( dfsIND[[i]], modelsIND[[i]], colIND, showZero=FALSE, showLine=TRUE)
	}

	legend( "topright", c('UnInduced', 'Induced  '), col=c(colUN,colIND), lwd=5, bg='white', cex=1.1)
	if ( ! is.null(showGrowth)) {
		pValText <- formatC( showPvalue, format="e", digits=2)
		legend( 'bottom', paste( "Delta Growth = ", round(showGrowth,digits=4), 
				"    P-value = ", pValText), bg='white', cex=1.1)
	}
	
	if (asPNG) {
		pngfile <- file.path( plot.path, paste( fileRoot, "png", sep="."))
		dev.print( png, pngfile, width=900, height=720)
	}
}


plotOneGene.TFKO <- function( gene, dfs, models, 
			plot.path=".", plot.prefix="", asPNG=FALSE, 
			showGrowth=NULL, showPvalue=NULL) {

	checkX11( bg="white", width=10, height=7)
	par( mai=c( 1,1,0.8,0.4))
	fileRoot <- paste( gene, plot.prefix, sep="_")

	# scan all the data frames for plot extents
	allDays <- allValues <- allClass <- vector()
	nModels <- length( dfs)
	if (nModels < 1) return()
	for (i in 1:nModels) {
		sml <- dfs[[i]]
		if ( is.null(sml)) next
		# small chance of extra non-drawn data to prevent breaking the model fitting
		if ( "DrawIt" %in% colnames(sml)) sml <- subset( sml, DrawIt)
		allDays <- c( allDays, sml$Day)
		allValues <- c( allValues, sml$Value)
	}
	xlim <- range( allDays, na.rm=TRUE) + c( -0.5,0.5)
	ylim <- range( allValues, na.rm=TRUE)
	ylim[2] <- ylim[2] + min( diff(ylim), 3) *0.35
	if ( ylim[1] > 0) ylim[1] <- ylim[1] * 0.9

	mainText <- paste( "Mutant TFKO Timecourse:    ", plot.prefix, "\n", gene, " -  ", gene2Product(gene))
	plot( 1, 1, type='n', main=mainText, xlim=xlim, ylim=ylim, xlab="Time  (days)", xaxt='n',
			ylab="Abundance   Log2(RPMHEG)", font.axis=2, font.lab=2)
	axis( side=1, at=sort( unique( allDays)), font=2)


	# local function to draw one "Linear Model" object...
	showOneTimeCourse <- function( df, lmAns, thisCol, thisPCH=21, pt.cex=1, showZero=TRUE, showLine=TRUE) {

		if ( is.null(df)) return()

		thisX <- df$Day
		thisY <- df$Value
		xFac <- factor( thisX)
		yMeans <- rep.int( 0, nlevels(xFac))
		theseCol <- rep.int( thisCol, length(thisX))
		
		# perhaps do not redraw the zero points (where 'Zero' includes day 1)
		toShow <- 1:length(thisX)
		if ( !showZero) toShow <- which( ! (thisX %in% DAY_ZERO))
		theseCol[ thisX %in% DAY_ZERO] <- 'black'

		points( jitter(thisX,factor=0.05)[toShow], thisY[toShow], pch=thisPCH, bg=theseCol[toShow], cex=pt.cex)
		yMeans <- tapply( thisY, xFac, mean)
		mX <- as.numeric( names(yMeans))
		for ( who in which( yMeans > 0)) {
			if ( !showZero && mX[who] %in% DAY_ZERO) next
			lines( mX[who]+c(-0.2,0.2), rep.int(yMeans[who],2), col=thisCol, lwd=2)
		}

		if ( is.null(lmAns)) return()
		if ( showLine && length(yMeans) > 1) {
			abline.segment( reg=lmAns, x1=min(mX), x2=max(mX), col=thisCol, lty=2, lwd=2)
		}
	}


	# OK, all set up, draw those time courses 
	colSet <- rainbow( nModels, end=0.72)
	col <- "seagreen"
	
	# only show the day zero one time, and always show the best fit 
	for ( i in 1:nModels) {
		if ( is.null( dfs[[i]])) next
		showZero <- (i == nModels)
		showOneTimeCourse( dfs[[i]], models[[i]], colSet[i], pt.cex=1.25, showZero=showZero, showLine=TRUE)
	}

	legend( "topright", names(dfs), col=colSet, lwd=3, bg='white', cex=1-(nModels*0.02))
	if ( ! is.null(showGrowth)) {
		pValText <- formatC( showPvalue, format="e", digits=2)
		legend( 'bottom', paste( "Delta Growth = ", round(showGrowth,digits=4), 
				"    P-value = ", pValText), bg='white', cex=1.1)
	}
	
	if (asPNG) {
		pngfile <- file.path( plot.path, paste( fileRoot, "png", sep="."))
		dev.print( png, pngfile, width=900, height=720)
	}
}


# model the growth defects and make plots and growth defect summarys for one sample key file
model.OneTFOEsample <- function( sampleKeyFile="SamplesKey.Run1.47TFlogPre.txt",
				sharedDayZeroUninduced=NULL, makePlots=TRUE) {

	# load that sample key
	allSamples <- loadSampleKeyFile( sampleKeyFile)
	nSamples <- nrow( allSamples)

	# the 'Run Name' is extracted from the name of the Samples Key filename for the results folder name
	results.path <- createResultsFolder( sampleKeyFile)

	# let's allow more than one 'Experiment' per Sample Key file
	expFactor <- factor( allSamples$Experiment)
	tapply( 1:nrow(allSamples), expFactor, function(x) {

		# do the subset of samples for this experiment
		samples <- allSamples[x, ]
		nSamples <- nrow( samples)

		# get the final results table for experiment
		myExperiment <- samples$Experiment[1]
		cat( "\nModelling TFOE growth defects for experiment: ", myExperiment)

		resultsFileIn <- paste( "FinalResults", myExperiment, "txt", sep=".")
		resultsFileIn <- file.path( results.path, resultsFileIn)
		if ( ! file.exists( resultsFileIn)) {
			cat( "\nResults file not found:  ", resultsFileIn)
			cat( "\nSkipping...")
			return()
		}
		tbl <- read.delim( resultsFileIn, as.is=T)
		tbl <- cleanMetaDataFields( tbl)

		# create folder for growth defect plots
		plot.path <- file.path( results.path, paste( "TFOE.GrowthDecayPlots", myExperiment, sep="."))
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=TRUE)

		# default behavior is to process each sub class separately, and merge them all back together
		out <- data.frame()
		facClass <- factor(tbl$Class)
		facSub1 <- factor(tbl$SubClass1)
		facSub2 <- factor(tbl$SubClass2)
		
		# allow a common shared set of Day Zero uninduced measrurements to be shared
		commonDayZeroUNDrows <- vector()
		if ( ! is.null( sharedDayZeroUninduced)) {
			# first place to check is SubClass1
			if ( any (tbl$SubClass1 == sharedDayZeroUninduced)) {
				commonDayZeroUNDrows <- which( tbl$Day %in% DAY_ZERO & tbl$SubClass1 == sharedDayZeroUninduced)
			}
			#if ( any (tbl$Class == sharedDayZeroUninduced)) {
			#	commonDayZeroUNDrows <- which( tbl$Day %in% DAY_ZERO & tbl$Class == sharedDayZeroUninduced)
			#}
		}
		
		by( tbl, INDICES=list( facClass, facSub1, facSub2), function(x) {

			# make sure we have data, and get the class details
			if (nrow(x) < 1) return()
			desc <- createClassDescriptor(x)
			cat( "\nModelling:  ", desc, "\n")
			myclass <- x$Class[1]
			mysub1 <- x$SubClass1[1]
			mysub2 <- x$SubClass2[1]
			
			# allow the addition of the common shared day zero uninduced data
			if ( !is.null(sharedDayZeroUninduced) && mysub1 != sharedDayZeroUninduced) {
				if ( length(commonDayZeroUNDrows)) {
					smlDF <- tbl[ commonDayZeroUNDrows, ]
					x <- rbind( x, smlDF)
				}
			}

			# there may be cases where all the rows in a SubClass are day 0/1. as when its a negative control
			# group.   If so, catch it here and do nothing
			if ( all (x$Day %in% DAY_ZERO)) return()

			# call the modeler
			ans <- model.TFOE.GrowthDefects(x, makePlots=makePlots, plot.path=plot.path,
						plot.prefix=desc, asPNG=TRUE)
			saveModelAns <<- ans

			# we get back matrices of Rate and Pvalues, linearize them
			rates <- ans$GrowthDefect
			pvals <- ans$P_Value
			days <- as.numeric( sub( "Day", "", rownames(rates)))
			genes <- colnames(rates)
			nDays <- length(days)
			nGenes <- length(genes)
			smlday <- smlgene <- smlrate <- smlpval <- vector()
			nout <- 0
			for ( i in 1:nGenes) {
				nnow <- (nout+1) : (nout+nDays)
				smlday[nnow] <- days
				smlgene[nnow] <- genes[i]
				smlrate[nnow] <- rates[,i]
				smlpval[nnow] <- pvals[,i]
				nout <- nout + nDays
			}
			plotID <- paste( smlgene, desc, sep="_")
			sml <- data.frame( "PlotID"=plotID, "Gene"=smlgene, "Product"=gene2Product(smlgene), 
					"Class"=myclass, "SubClass1"=mysub1,
					"SubClass2"=mysub2, "Day"=smlday, "Delta_Growth"=round(smlrate,digits=4),
					"Delta_Pvalue"=smlpval, stringsAsFactors=FALSE)

			out <<- rbind( out, sml)
		})

		# when all done, we can sort and write out
		ord <- diffExpressRankOrder( -out$Delta_Growth, out$Delta_Pvalue)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)

		resultsFileOut <- paste( "TFOE.GrowthDefect", myExperiment, "txt", sep=".")
		resultsFileOut <- file.path( results.path, resultsFileOut)
		write.table( out, resultsFileOut, sep="\t", quote=F, row.names=F)

		# for the HTML, drop meta data columns with no data
		if ( all( out$SubClass2 == "")) out <- out[ , -match("SubClass2",colnames(out))]
		if ( all( out$SubClass1 == "")) out <- out[ , -match("SubClass1",colnames(out))]
		htmlFile <- sub( "txt$", "html", resultsFileOut)
		colnames(out) <- sub( "_", " ", colnames(out))
		table2html( out, htmlFile, title=paste( "TFOE Growth Defect Results:   ", myExperiment),
			linkColumnNames="PlotID", linkPaths=basename(plot.path))

	})  ## done with tapply() of all experiments in the run

	cat( "\nDone.\n")
}


# model the growth defects and make plots and growth defect summarys for one sample key file
model.OneTFKOsample <- function( sampleKeyFile="SamplesKey.Run8.TFKO.txt", makePlots=TRUE) {

	# load that sample key
	allSamples <- loadSampleKeyFile( sampleKeyFile)
	nSamples <- nrow( allSamples)

	# the 'Run Name' is extracted from the name of the Samples Key filename for the results folder name
	results.path <- createResultsFolder( sampleKeyFile)

	# let's allow more than one 'Experiment' per Sample Key file
	expFactor <- factor( allSamples$Experiment)
	tapply( 1:nrow(allSamples), expFactor, function(x) {

		# do the subset of samples for this experiment
		samples <- allSamples[x, ]
		nSamples <- nrow( samples)

		# get the final results table for experiment
		myExperiment <- samples$Experiment[1]
		cat( "\nModelling TFKO growth defects for experiment: ", myExperiment)

		resultsFileIn <- paste( "FinalResults", myExperiment, "txt", sep=".")
		resultsFileIn <- file.path( results.path, resultsFileIn)
		if ( ! file.exists( resultsFileIn)) {
			cat( "\nResults file not found:  ", resultsFileIn)
			cat( "\nSkipping...")
			return()
		}
		tbl <- read.delim( resultsFileIn, as.is=T)
		tbl <- cleanMetaDataFields( tbl)

		# create folder for growth defect plots
		plot.path <- file.path( results.path, paste( "TFKO.GrowthDecayPlots", myExperiment, sep="."))
		if ( ! file.exists( plot.path)) dir.create( plot.path, recursive=TRUE)

		# default behavior is to process each sub class separately, and merge them all back together
		out <- data.frame()
		facClass <- factor(tbl$Class)
		facSub1 <- factor(tbl$SubClass1)
		facSub2 <- factor(tbl$SubClass2)
		
		# allow a common shared set of Day Zero uninduced measrurements to be shared
		# rather, for TFKO, all models alwasy include Day0
		commonDayZeroUNDrows <- which( tbl$Day %in% DAY_ZERO)


		by( tbl, INDICES=list( facClass, facSub1, facSub2), function(x) {

			# make sure we have data, and get the class details
			if (nrow(x) < 1) return()
			desc <- createClassDescriptor(x)
			cat( "\nModelling:  ", desc, "\n")
			myclass <- x$Class[1]
			mysub1 <- x$SubClass1[1]
			mysub2 <- x$SubClass2[1]
			
			# force the addition of the common shared day zero uninduced data
			# already done at the calling routine -- don't add them in again
			#smlDF <- tbl[ commonDayZeroUNDrows, ]
			#x <- rbind( x, smlDF)

			# there may be cases where all the rows in a SubClass are day 0/1. as when its a negative control
			# group.   If so, catch it here and do nothing
			if ( all (x$Day %in% DAY_ZERO)) return()

			# call the modeler
			ans <- model.TFKO.GrowthDefects(x, makePlots=makePlots, plot.path=plot.path,
						plot.prefix=desc, asPNG=TRUE)
			saveModelAns <<- ans

			# we get back matrices of Rate and Pvalues, linearize them
			rates <- ans$GrowthDefect
			pvals <- ans$P_Value
			days <- as.numeric( sub( "Day", "", rownames(rates)))
			genes <- colnames(rates)
			nDays <- length(days)
			nGenes <- length(genes)
			smlday <- smlgene <- smlrate <- smlpval <- vector()
			nout <- 0
			for ( i in 1:nGenes) {
				nnow <- (nout+1) : (nout+nDays)
				smlday[nnow] <- days
				smlgene[nnow] <- genes[i]
				smlrate[nnow] <- rates[,i]
				smlpval[nnow] <- pvals[,i]
				nout <- nout + nDays
			}
			plotID <- paste( smlgene, desc, sep="_")
			sml <- data.frame( "PlotID"=plotID, "Gene"=smlgene, "Product"=gene2Product(smlgene), 
					"Class"=myclass, "SubClass1"=mysub1,
					"SubClass2"=mysub2, "Day"=smlday, "Delta_Growth"=round(smlrate,digits=4),
					"Delta_Pvalue"=smlpval, stringsAsFactors=FALSE)

			out <<- rbind( out, sml)
		})

		# when all done, we can sort and write out
		ord <- diffExpressRankOrder( -out$Delta_Growth, out$Delta_Pvalue)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)

		resultsFileOut <- paste( "TFKO.GrowthDefect", myExperiment, "txt", sep=".")
		resultsFileOut <- file.path( results.path, resultsFileOut)
		write.table( out, resultsFileOut, sep="\t", quote=F, row.names=F)

		# for the HTML, drop meta data columns with no data
		if ( all( out$SubClass2 == "")) out <- out[ , -match("SubClass2",colnames(out))]
		if ( all( out$SubClass1 == "")) out <- out[ , -match("SubClass1",colnames(out))]
		htmlFile <- sub( "txt$", "html", resultsFileOut)
		colnames(out) <- sub( "_", " ", colnames(out))
		table2html( out, htmlFile, title=paste( "TFKO Growth Defect Results:   ", myExperiment),
				linkColumnNames="PlotID", linkPaths=basename(plot.path))

	})  ## done with tapply() of all experiments in the run

	cat( "\nDone.\n")
}


poolHitsSummary <- function( folder, bowtieMode=c("Bowtie2","Bowtie1")) {

	# generate a summary of valid read counts that hit the expected pool of TFI genes
	# versus false hits to other genes not in the pool
	# 'folder' is one results folder from one run of 'align.OneSample()'
	bowtieMode <- match.arg( bowtieMode)

	allPattern <- paste( "^AllGenes.RawCounts.*", bowtieMode, "txt$", sep=".")
	allFiles <- dir( folder, pattern=allPattern, full.name=T)
	poolFiles <- sub( "/All", "/Pool", allFiles)

	experiments <- sub( "(AllGenes.RawCounts.)(.+)(.Bowtie[12].txt)", "\\2", basename(allFiles))
	nExp <- length( experiments)

	allCounts <- poolCounts <- expNames <- colNames <- vector()
	for ( i in 1:length( allFiles)) {
		allTbl <- read.delim( allFiles[i], as.is=T)
		poolTbl <- read.delim( poolFiles[i], as.is=T)
		if ( ncol(allTbl) != ncol(poolTbl)) {
			cat( "\nError:  column count mismatch of 'All' & 'Pool' results: ", allFiles[i])
			next
		}
		if ( any( colnames(allTbl) != colnames(poolTbl))) {
			cat( "\nError:  column names mismatch of 'All' & 'Pool' results: ", allFiles[i])
			next
		}

		# get the column sums for each, and accumulate
		myCntsAll <- apply( as.matrix( allTbl[ ,2:ncol(allTbl)]), 2, sum, na.rm=T)
		myCntsPool <- apply( as.matrix( poolTbl[ ,2:ncol(poolTbl)]), 2, sum, na.rm=T)
		allCounts <- c( allCounts, myCntsAll)
		poolCounts <- c( poolCounts, myCntsPool)
		expNames <- c( expNames, rep.int( experiments[i], length(myCntsAll)))
		colNames <- c( colNames, names(myCntsAll))
	}

	# now we can know the percent that hit the pool
	meanAll <- round( mean( allCounts, na.rm=T))
	meanPool <- round( mean( poolCounts, na.rm=T))
	pctPool <- poolCounts * 100 / allCounts
	meanPct <- round( mean( pctPool, na.rm=T), digits=2)
	sdPct <- round( sd( pctPool, na.rm=T), digits=2)

	# report results to stdout
	cat( "\nFolder:      ", folder)
	cat( "\nExperiments: ", experiments)
	cat( "\n\nOverall Pool Hit Summary:")
	cat( "\nMean 'All'  Read Count per Column:  ", prettyNum(meanAll,big.mark=','))
	cat( "\nMean 'Pool' Read Count per Column:  ", prettyNum(meanPool,big.mark=','))
	cat( "\nMean % of Reads hitting Pool Genes: ", meanPct)
	cat( "\nStandard Deviation of % Pool:       ", sdPct)
	cat( "\n\nPercentage by Experiment:\n")
	print( tapply( pctPool, factor(expNames), function(x) return( round( mean( x, na.rm=T), digits=2))))
}
