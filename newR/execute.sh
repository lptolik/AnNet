#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1/datasets

##---algorithms available in Rclustering.R
#alg     <- vector(length=16)
#alg[1]  <- "louvain"
#alg[2]  <- "infomap"
#alg[3]  <- "fc"
#alg[4]  <- "lec"
#alg[5]  <- "sgG1"
#alg[6]  <- "sgG2"
#alg[7]  <- "sgG5"
#alg[8]  <- "wt"
#alg[9]  <- "spectral"
#alg[10] <- "louvain2"
#alg[11] <- "infomap2"
#alg[12] <- "fc2"
#alg[13] <- "lec2"
#alg[14] <- "sgG12"
#alg[15] <- "wt2"
#alg[16] <- "Spectral2"
##---

## parameters for subsampling.R script
ALG=10
TYPE=2
MASK=20
CNMIN=-1
CNMAX=-1

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

## load subsampling.R script
cp -r $SCRIPTDIR/Rclustering.R .

## load graph
cp -r $DATADIR/PPI_Presynaptic.gml .

## initiallise environment module
. /etc/profile.d/modules.sh

## load module R
module load R 

## run subsampling on clustering algorithm
time Rscript Rclustering.R $ALG $TYPE $MASK $CNMIN $CNMAX

## define output loocation for subsampling consensus file
OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.txt $OUTDIR

echo "$0 done!"

exit 0
