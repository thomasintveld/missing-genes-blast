missing-genes-blast
===================

Manual operations to perform at this point:

* Install a command-line version of the ncbi blast tools. Especially 'blastp' needs to run from the commandline in any location.
* Install the protein databases in blast for homo, gallus, taeniopygia, xenopus, anolis and latimeria. These need to be run through the create-database script in blast.

The relevant protein FASTA files can be downloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html):
* [Homo Sapiens](http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz) (release 75)
* [Gallus Gallus](http://ftp.ensembl.org/pub/release-75/fasta/gallus_gallus/pep/Gallus_gallus.Galgal4.75.pep.all.fa.gz) (chicken, release 75)
* [taeniopygia](http://ftp.ensembl.org/pub/release-75/fasta/taeniopygia_guttata/pep/Taeniopygia_guttata.taeGut3.2.4.75.pep.all.fa.gz) (zebrafinch, release 75)
* [anolis carolinensis](http://ftp.ensembl.org/pub/release-75/fasta/anolis_carolinensis/pep/Anolis_carolinensis.AnoCar2.0.75.pep.all.fa.gz) (lizard, release 75)
* [xenopus tropicalis](http://ftp.ensembl.org/pub/release-75/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.JGI_4.2.75.pep.all.fa.gz) (frog, release 75)
* [latimeria](http://ftp.ensembl.org/pub/release-75/fasta/latimeria_chalumnae/pep/Latimeria_chalumnae.LatCha1.75.pep.all.fa.gz) (coelacanth fish, release 75)

Get these needed FASTA files from ncbi (protein files) for the six species, put them in the data/fasta folder and make sure the locations are linked correctly in the autoBLAST.R scripts. 

Apart from that, everything should be automised: 
* First run the 10 (5 species both ways) autoBLAST.R scripts. These generate full lists of found/not found proteins in the reciprocal species, and put them in the data/results_blast folder.
* Then run the twowayAnalysis.R script. This runs an interpretation layer on top of the BLAST results, to perform a proper two-way reciprocal analysis. This automatically merges all the results together to get a final merged file.

Author:
======
Thomas in't Veld (KU Leuven, November 2014).

Contributors:
=============
- Frans Schuit
- Stefanie De Coster
- Lieven Thorrez

