
library(XML)

URLPRE <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="

load('../data/genelists/infohomo.RData')
infohomo <- infosubject

genes.with.names <- data.frame(ensgid=infohomo$Ensembl.Gene.ID, assname=infohomo$Associated.Gene.Name, chromosome=infohomo$Chromosome.Name)
genes.with.names <- unique(genes.with.names)

# limit to our interest in chromosomes 7 and 17
genes.with.names <- genes.with.names[genes.with.names$chromosome %in% c(7, 17),]

number.pubmed.checker <- function(gene.name){
  
  url <- paste(URLPRE, gene.name, sep='')
  
  xml <- xmlParse(url)
  data <- xmlToList(xml)
  if(length(data)>0){
    count <- as.numeric(as.character(data$Count))
  } else {
    count <- 0
  }	
  
  return(count)
  
}

df <- sapply(genes.with.names$assname, function(x) number.pubmed.checker(x))
result <- data.frame(ass.name=genes.with.names$assname, pub.count=as.numeric(df))
write.csv(result, '../publication_count_chr7_17.csv', row.names=FALSE)
