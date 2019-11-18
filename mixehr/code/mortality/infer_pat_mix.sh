#!/bin/bash

# train mixehr on the full mimic dataset for 10000 iterations or until convergence within 10^-6
k=75
mixehr=/code/mixehr/mixehr

outdir=/results/mortality
if [ ! -d $outdir ]; then mkdir -p $outdir; fi

data0=/data/Mimic/version_1_3/processed/mimic_trainData.txt
meta0=/data/Mimic/version_1_3/processed/mimic_meta.txt
testdata0=/data/Mimic/version_1_3/processed/mimic_testData_firstAdm.txt
data=$outdir/$(basename $data0)
meta=$outdir/$(basename $meta0)
testdata=$outdir/$(basename $testdata0)
ln - s $data0 $data
ln - s $meta0 $meta
ln - s $testdata0 $testdata

testdata=$outdir/testData_fold${cvfold}.txt

niter=500

# train mixehr on training fold
$mixehr -f $data -m $meta -i $niter -k $k -n JCVB0

# get the iter in case converged before $niter
trainedPrefix=$(ls /results/mimic_trainData_JCVB0_nmar_K${k}_iter*_phi.csv | sed s/_phi.csv//)

# infer metaphe on both training and testing fold
$mixehr -m $meta -n JCVB0 --newPatsData $data --trainedModelPrefix $trainedPrefix -k $k --inferNewPatentMetaphe --inferPatParams_maxiter 100
$mixehr -m $meta -n JCVB0 --newPatsData $testdata --trainedModelPrefix $trainedPrefix -k $k --inferNewPatentMetaphe --inferPatParams_maxiter 100