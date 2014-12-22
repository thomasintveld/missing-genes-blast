library(ggplot2)

load.frame <- function(file){
	load(file)
	return(result)
	#return(read.csv(file, stringsAsFactors=FALSE))

}

append.and.order.by.location <- function(df){
	
	load('InfoPerGene.RData')
  info <- infopergene

	#info <- read.csv('InfoPerGene.csv', stringsAsFactors=FALSE)

	all <- merge(df, info, by='Ensembl.Gene.ID')
	
	all[all$Chromosome.Name=='1',]['Chromosome.Name'] <- '01'
	all[all$Chromosome.Name=='2',]['Chromosome.Name'] <- '02'
	all[all$Chromosome.Name=='3',]['Chromosome.Name'] <- '03'
	all[all$Chromosome.Name=='4',]['Chromosome.Name'] <- '04'
	all[all$Chromosome.Name=='5',]['Chromosome.Name'] <- '05'
	all[all$Chromosome.Name=='6',]['Chromosome.Name'] <- '06'
	all[all$Chromosome.Name=='7',]['Chromosome.Name'] <- '07'
	all[all$Chromosome.Name=='8',]['Chromosome.Name'] <- '08'
	all[all$Chromosome.Name=='9',]['Chromosome.Name'] <- '09'
	

	all <- all[order(all$Chromosome.Name, all$Gene.Start..bp.),]

	return(all)

}

create.rolling.col <- function(df, window.size){

	# missing column is $missing

	windows <- 1:length(df$missing)

	rolling.sum <- sapply(windows, function(i) if(i >= window.size){ return(sum(df$missing[(i - window.size + 1):i])) }else{ return(0)} )

	df['rolling.sum'] <- rolling.sum
	df['index'] <- 1:nrow(df)

	return(df)

}

plot.rolling <- function(df){

	plot <- ggplot(data=df) + geom_line(aes(x=index, y=rolling.sum), color='maroon4', size=0.7) + ylab('Rolling Window Sum') + xlab('Location Index')
	
	return(plot)

}

do.all <- function(file, window.size){

	df <- load.frame(file)
	df <- append.and.order.by.location(df)
	df <- create.rolling.col(df, window.size)
	df <- merge.method(df)
	
	return(df)


}

merge.with.manual <- function(file, window.size){
	
	df <- do.all(file, window.size)
	
	#df <- merge.method(df)
	
	return(df)
	
}

plot.file <- function(file, window.size){

	df <- load.frame(file)
	df <- append.and.order.by.location(df)
	df <- create.rolling.col(df, window.size)

	return(plot.rolling(df))


}

mark.peaks <- function(df, peak.threshold){
	
	df['peak.auto'] <- sapply(df$rolling.sum, function(rolsum) if(is.null(rolsum)){return(NA)} else if(rolsum>peak.threshold){return('auto')}else{return(NA)})
	
	if(!is.null(df$Manual.Sum)){
		df['peak.manual'] <- sapply(df$Manual.Sum, function(rolsum) 
			if(!is.na(rolsum) & rolsum>peak.threshold){return('manual')}else{return(NA)} )
	
	}

	return(df)
}

plot.peaks <- function(df){
	
	plotbase <- ggplot(data=df)
	p1 <- geom_line(aes(x=index, y=peak.auto), color='maroon', size=50)
	p2 <- geom_line(aes(x=index, y=peak.manual), color='darkslategrey', size=50)
	
	return(plotbase + p1 + p2 + scale_y_continuous(limit=c(0.95, 1.15)))
	
}

plot.with.manual <- function(df){
	plotbase <- ggplot(data=df)
	p1 <- geom_line(aes(x=index, y=rolling.sum), color='maroon4', size=0.6)
	p2 <- geom_line(aes(x=index, y=Manual.Sum), color='darkslategrey')
	
	return(plotbase + p1 + p2 + ylab('Rolling Window Sum') + xlab('Location Index'))

}

merge.method <- function(df){
	
	manual <- read.csv('../../manuele_missinggenes_sums.csv', stringsAsFactors=FALSE)
	
	all <- merge(df, manual, by='Associated.Gene.Name', all.x=TRUE, all.y=FALSE)
	all <- all[order(all$Chromosome.Name, all$Gene.Start..bp.),]
	
	return(all)

}

compare.with.manual.plot <- function(df){


	p1 <- plot.rolling(df)
	
	p2 <- ggplot(data=df) + geom_line(aes(x=index, y=Manual.Sum), color='darkslategrey', size=0.7) + ylab('Rolling Window Sum (Manual)') + xlab('Location Index')
	
	
	
	pdf('../../20140914 - presentation thomas/presentation_threshold50peaks_withmanual.pdf', width=14, height=10)
	multiplot(p2, p1, cols=1)
	dev.off()


}

multiple.thresholds <- function(window){

	p1 <- do.all('resultTwoWayThreshold652014-07-17 17:14:09.csv', window) + ggtitle(paste('Missing Threshold 65% (rolling.window ', window,')', sep=''))
	p2 <- do.all('resultTwoWayThreshold702014-07-11 18:48:07.csv', window) + ggtitle(paste('Missing Threshold 70% (rolling.window ', window,')', sep=''))
	p3 <- do.all('resultTwoWayThreshold752014-06-28 04:32:30.csv', window) + ggtitle(paste('Missing Threshold 75% (rolling.window ', window,')', sep=''))

	pdf(paste("ThreeThresholdsWindow", window, ".pdf", sep=''), width=7, height=10)
	multiplot(p1, p2, p3, cols=1)
	dev.off()


}

multiple.thresholds.file <- function(file){

  p1 <- do.all(file, 50) + ggtitle(paste('Missing Threshold 65% (rolling.window 50)', sep=''))
  p2 <- do.all(file, 100) + ggtitle(paste('Missing Threshold 65% (rolling.window 100)', sep=''))
  p3 <- do.all(file, 200) + ggtitle(paste('Missing Threshold 65% (rolling.window 200)', sep=''))

  pdf(paste("../../ThreeThresholdsWindow65.pdf", sep=''), width=7, height=10)
  multiplot(p1, p2, p3, cols=1)
  dev.off()


}

graph.peaks.thresholds <- function() {
	
	df35 <- mark.peaks(merge.with.manual('resultTwoWayThreshold352014-09-06.csv', 200), 15)
	df40 <- mark.peaks(merge.with.manual('resultTwoWayThreshold402014-09-05.csv', 200), 15)
	df45 <- mark.peaks(merge.with.manual('resultTwoWayThreshold452014-09-04.csv', 200), 15)
	df50 <- mark.peaks(merge.with.manual('resultTwoWayThreshold502014-09-03.csv', 200), 15)
	df55 <- mark.peaks(merge.with.manual('resultTwoWayThreshold552014-09-03.csv', 200), 15)
	df60 <- mark.peaks(merge.with.manual('resultTwoWayThreshold602014-08-29.csv', 200), 15)
	df63 <- mark.peaks(merge.with.manual('resultTwoWayThreshold632014-09-01.csv', 200), 15)
	df65 <- mark.peaks(merge.with.manual('resultTwoWayThreshold652014-08-29 05:22:26.csv', 200), 15)
	df70 <- mark.peaks(merge.with.manual('resultTwoWayThreshold702014-08-29.csv', 200), 15)
	
	
	peak.manual <- data.frame(
		index=df60$index,
		chromosome=df60$Chromosome.Name,
		start.bp=df60$Gene.Start..bp.,
		peak=df60$peak.manual,
		type=as.character('manual')
	)
	
	peak.35 <- data.frame(
		index=df35$index,
		chromosome=df35$Chromosome.Name,
		start.bp=df35$Gene.Start..bp.,
		peak=df35$peak.auto,
		type=as.character(35)
	)
	
	peak.40 <- data.frame(
		index=df40$index,
		chromosome=df40$Chromosome.Name,
		start.bp=df40$Gene.Start..bp.,
		peak=df40$peak.auto,
		type=as.character(40)
	)
	
	peak.45 <- data.frame(
		index=df45$index,
		chromosome=df45$Chromosome.Name,
		start.bp=df45$Gene.Start..bp.,
		peak=df45$peak.auto,
		type=as.character(45)
	)
	
	peak.50 <- data.frame(
		index=df50$index,
		chromosome=df50$Chromosome.Name,
		start.bp=df50$Gene.Start..bp.,
		peak=df50$peak.auto,
		type=as.character(50)
	)
	
	peak.55 <- data.frame(
		index=df55$index,
		chromosome=df55$Chromosome.Name,
		start.bp=df55$Gene.Start..bp.,
		peak=df55$peak.auto,
		type=as.character(55)
	)
	
	peak.60 <- data.frame(
		index=df60$index,
		chromosome=df60$Chromosome.Name,
		start.bp=df60$Gene.Start..bp.,
		peak=df60$peak.auto,
		type=as.character(60)
	)
	
	peak.63 <- data.frame(
		index=df63$index,
		chromosome=df63$Chromosome.Name,
		start.bp=df63$Gene.Start..bp.,
		peak=df63$peak.auto,
		type=as.character(63)
	)
	
	peak.65 <- data.frame(
		index=df65$index,
		chromosome=df65$Chromosome.Name,
		start.bp=df65$Gene.Start..bp.,
		peak=df65$peak.auto,
		type=as.character(65)
	)
	
	peak.70 <- data.frame(
		index=df70$index,
		chromosome=df70$Chromosome.Name,
		start.bp=df70$Gene.Start..bp.,
		peak=df70$peak.auto,
		type=as.character(70)
	)
	
	all.peaks <- rbind(peak.manual, peak.50, peak.60, peak.65, peak.70)
	#all.peaks <- all.peaks[!is.na(all.peaks$peak),]
	
	p <- ggplot(data=all.peaks, aes(x=index, xend=index+8, y=chromosome, yend=chromosome	, colour=factor(peak))) + geom_segment(size=20) + ylab('Threshold') + xlab('chrom') + guides(color=FALSE) + scale_color_manual(values=c('darkslategrey', 'maroon4')) + scale_x_discrete(aes(x=factor(chromosome)))
	
	q <- ggplot(data=all.peaks, aes(x=index, xend=index+1, y=type, yend=type, colour=factor(chromosome))) + geom_segment(size=20) + ylab('Threshold') + xlab('Location') + guides(color=FALSE) + scale_color_manual(values=c('darkslategrey', 'maroon4'))			
				
		p <- ggplot(data=all.peaks, aes(x=index, xend=index+8, y=chromosome, yend=chromosome	, colour=factor(peak))) + geom_segment(size=10) + ylab('Chromosome') + xlab('index') + guides(color=FALSE) + scale_color_manual(values=c('darkslategrey', 'maroon4')) + scale_x_discrete(aes(x=factor(chromosome)))			
				
	return(p)
}

banding.one.chromosome <- function(chromosome){
	df50 <- mark.peaks(merge.with.manual('resultTwoWayThreshold502014-09-03.csv', 200), 15)
	
	peak.50 <- data.frame(
		index=df50$index,
		chromosome=df50$Chromosome.Name,
		start.bp=df50$Gene.Start..bp.,
		peak=df50$peak.auto,
		type=as.character(50)
	)
	
	all.peaks <- peak.50[peak.50$chromosome==chromosome,]
	
	p <- ggplot(data=all.peaks, aes(x=index, xend=index+3, y=chromosome, yend=chromosome	, colour=factor(peak))) + geom_segment(size=10) + ylab('Chromosome') + xlab('index')+ guides(color=FALSE)
	
	return(p)
	
	

}



#### config


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

analysePeaks <- function(df, threshold){
	
	# fix rownames of df
	rownames(df) <- 1:nrow(df)
	
	
	df <- df[df$rolling.sum>threshold,]
	indices <- 1:nrow(df)
	
		
	delta <- sapply(indices, function(i) if(i==1){return(1)}else{return( as.integer(rownames(df[i,])) - as.integer(rownames(df[i-1,]))  )})
	df['deltas'] <- delta
	df[df$deltas>1,]['deltas'] <- 2
	df['deltas'] <- df$deltas - 1
	peak <- sapply(indices, function(i) return( sum(df$deltas[1:i])+1))
	df['peak'] <- peak
	#df['index'] <- rownames(df)
	
	df <- split(df, df$peak)
	
	average.position <- sapply(df, function(peak) return(mean(peak$Gene.Start..bp.)))
	chromosomes <- sapply(df, function(peak) return(table(peak$Chromosome.Name)))
	average.exon.gc.content <- sapply(df, function(peak) return(mean(peak$GC.exon.count)))
	average.length <- sapply(df, function(peak) return(mean(peak$length)))
	number.of.genes <- sapply(df, function(peak) return(nrow(peak)))
	average.index <- sapply(df, function(peak) return(mean(peak$index)))
	peak <- 1:length(df)
	
	infotable <- cbind(peak, average.position, average.exon.gc.content, average.length, number.of.genes)
	
	return(infotable)

}

analyseGCcontent.naive <- function(df, N){

    # idea is to build a normal distribution out of the GC-averages of 1000 samples of 631 genes,
    # and then to compare the GC average of the missing genes with that distribution.
	nrOfMissing <- nrow(df[df$missing,])

	indices <- 1:N
	distribution <- sapply(indices, function(i) getAverageGC(df[sample(1:nrow(df), nrOfMissing),]))
	
	meanDist <- mean(distribution)
	sigmaDist <- sd(distribution)
	
	meanMissing <- getAverageGC(df[df$missing,])
	
	z <-  (meanMissing - meanDist)/(sigmaDist/sqrt(N))
	p <- 2*pnorm(-abs(z))
	
	return(p)

}

analyseGCContent.allStuff <- function(df,N){

	indices <- 1:N

	distribution <- sapply(indices, function(i) return(analyseGCContent.welchttest(df)))

	return(t(distribution))

}

analyseGCContent.welchttest <- function(df){

	nrOfMissing <- nrow(df[df$missing,])

	dfMissing <- df[df$missing,]
	dfAll <- df[sample(1:nrow(df), nrOfMissing),]

	test <- t.test(dfMissing$GC.exon.count, dfAll$GC.exon.count)

	return(c(pvalue=test$p.value, delta=abs(test$estimate[1] - test$estimate[2]))) # this will return the p-value and the delta between the two means.

}

getAverageGC <- function(df, weighted=TRUE){
	
	if(!weighted){
		return(mean(df$GC.exon.count))
	}
	
	df['weightedGC'] <- df$length*df$GC.exon.count
	
	average <- sum(df$weightedGC)/sum(df$length)
	
	return(average)
	

}

testingTheProblematicOnes <- function(df, problematic){

	results <- sapply(problematic, function(gene) return(df[df$Associated.Gene.Name==gene,]['missing']))
	results <- do.call('c', results)
	return(results)

}

testingTheProblematicOnesChicken <- function(df, problematic){

	results <- sapply(problematic, function(gene) return(df[df$Associated.Gene.Name==gene,]['gallus']))
	results <- do.call('c', results)
	return(results)

}

keeping.appearances <- function(df){
	
	index <- 1:nrow(df)
	
	sapply(index, function(i) print(paste('Testing', df[i,'Associated.Gene.Name'], '... gallus:', df[i,'gallus'], '... taeniopygia:', df[i,'taeniopygia'], '... latimeria:', df[i,'latimeria'], '... anolis:', df[i,'anolis'], '... xenopus:', df[i,'xenopus'], 'so MISSING: ', df[i,'missing'])))
data.frame(df$Associated.Gene.Name, df$gallus, df$taeniopygia, df$f$xenopus, df$anolis, df$missing)
}

gc.content.analysis <- function(df){

	set.missing <- df[df$missing,]
	set.peak <- set.missing[set.missing$rolling.sum>15,]
	
	sample.normal <- df[sample(1:nrow(df), 500),]
	sample.peak <- set.peak[sample(1:nrow(set.peak), 500),]
	
	sample.normal['marker'] <- 1
	sample.peak['marker'] <- 2
	
	all <- rbind(sample.normal, sample.peak)
	
	p <- ggplot(data=all, aes(x=GC.exon.count, fill=as.factor(marker)) ) + geom_histogram() 
	
	p <- ggplot(data=sample.normal) + geom_histogram(aes(x=GC.exon.count), fill='slategrey' ) + xlim(0.3, 0.8) + ylim(0, 80) 
			
	q <- ggplot(data=sample.peak) + geom_histogram(aes(x=GC.exon.count), fill='maroon4' ) + xlim(0.3, 0.8) + ylim(0, 80)
				
	
	p <- ggplot() + geom_histogram(aes(x=sample.normal$GC.exon.count), fill='slategrey' ) + geom_histogram(aes(x=sample.peak$GC.exon.count), fill='maroon4')
				+ geom_histogram(aes(data=sample.peak, x=GC.exon.count), color='maroon4')
	
	return(p)
	
}

length.analysis <- function(df){

	set.missing <- df[df$missing,]
	set.peak <- set.missing[set.missing$rolling.sum>15,]
	
	sample.normal <- df[sample(1:nrow(df), 500),]
	sample.peak <- set.peak[sample(1:nrow(set.peak), 500),]
	
	sample.normal['marker'] <- 1
	sample.peak['marker'] <- 2
	
	all <- rbind(sample.normal, sample.peak)
	
	p <- ggplot(data=all, aes(x=length, fill=as.factor(marker)) ) + geom_histogram() 
	
	p <- ggplot(data=sample.normal) + geom_histogram(aes(x=log(length)), fill='slategrey' ) + xlim(5, 15) + ylim(0, 70)
			
	q <- ggplot(data=sample.peak) + geom_histogram(aes(x=log(length)), fill='maroon4' ) + xlim(5,15) + ylim(0,70)
				
	
	p <- ggplot() + geom_histogram(aes(x=sample.normal$GC.exon.count), fill='slategrey' ) + geom_histogram(aes(x=sample.peak$GC.exon.count), fill='maroon4')
				+ geom_histogram(aes(data=sample.peak, x=GC.exon.count), color='maroon4')
	
	return(p)
	
}

