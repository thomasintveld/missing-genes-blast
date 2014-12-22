library(plyr)
library(dplyr)
source('interpretResultWithThreshold.R')

window.size <- 100
fileNameOld <- '../data/results_merged/result_45_2014-11-17.RData'
fileNameNew <- '../data/results_merged/result_50_2014-12-21.RData'

df <- load.frame(fileNameNew)
df <- append.and.order.by.location(df)
df <- create.rolling.col(df, window.size)

plot <- ggplot(data=df) + geom_line(aes(x=index, y=rolling.sum), color='maroon4', size=0.5) 
plot <- plot + ylab('Rolling Window Sum') + xlab('Location Index')
plot

plotOld <- plot
plotNew <- plot


chromosomes <- unique(df$Chromosome.Name)

sums <- sapply(chromosomes, function(chrx) sum(df[df$Chromosome.Name==chrx,'missing']))
sumsPerChromosome <- data.frame(chromosome=chromosomes, nrMissing=sums)

# keep only 'usual' chromosomes
usualChromosomes <- c(paste(0, 1:9, sep=''), 10:22, 'X', 'Y')

sumsFiltered <- sumsPerChromosome[sumsPerChromosome$chromosome %in% usualChromosomes,]

multiplot(plotOld, plotNew)
# plot this
plotSumsPerChromosomes2 <- ggplot(sumsFiltered) + geom_bar(aes(x=chromosome, y=nrMissing), stat='identity')

