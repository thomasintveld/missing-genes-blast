missing-genes-blast
===================

Manual operations to perform at this point:

* Install a command-line version of the ncbi blast tools. Especially 'blastp' needs to run from the commandline in any location.
* Install the protein databases in blast for homo, gallus, taeniopygia, xenopus, anolis and latimeria. These need to be run through the create-database script in blast.

(A fully installed version of blast with databases can be downloaded from https://www.dropbox.com/s/i9a5rkgwp7s5ovk/blast.tar.gz?dl=0 at this point. Just extract the folder in any location you want and add the path to $PATH).

* Get the needed FASTA files from ncbi (protein files) for the six species, put them in the data/fasta folder and make sure the locations are linked correctly in the autoBLAST.R scripts. 

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

