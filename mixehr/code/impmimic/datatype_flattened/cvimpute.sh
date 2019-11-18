#!/bin/bash

mixehr=/results/mixehr

k=75
niter=200

Rscript /code/impmimic/processCrossValidSets.R

outdir=/results/impmimic/datatype_flattened

if [ ! -d $outdir ]; then mkdir -p $outdir; fi

meta=$outdir/metainfo.txt

knn=100

for cvfold in {1..5}; do 

    data=$outdir/train${cvfold}.txt
    testdata=$outdir/test${cvfold}.txt

    # train on cvfold
    $mixehr -f $data -m $meta -i $niter -k $k --inferenceMethod JCVB0
    
    # infer patient mixture on this fold    
    $mixehr --metaFile $meta --topics $k \
            --trainDataFile $data \
            --inferTrainPatientMetaphe \
            --output_dir $outdir/impute${cvfold}\
            --trainPatMetapheFile train_pat_mix.csv\
            --trainedModelPrefix train${cvfold}_JCVB0_nmar_K${k}_iter${iter}
            
    # impute by knn
    $mixehr --metaFile $meta --topics $k \
            --trainDataFile $data \
            --imputeTargetsFile /results/impute_target_pheId.txt \
            --imputePatDataFile $testdata \
            --trainPatMetapheFile $outdir/impute${cvfold}/train_pat_mix.csv \
            --trainPatIdFile $outdir/impute${cvfold}/train_pat_mix_patId.csv \
            --knn_impute $knn \
            --output_dir $outdir/impute${i}_knn${knn} \
            --trainedModelPrefix train${i}_JCVB0_nmar_K${k}_iter${iter}
done
