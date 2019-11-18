#!/bin/bash

# compile mixehr software C++ code (require armadillo and boost)
# cd mixehr
# make

# train mixehr on the full mimic dataset for 10000 iterations or until convergence within 10^-6
# k=75
# mixehr=/code/mixehr/mixehr

# outdir=/results/mixmimic
# if [ ! -d $outdir ]; then mkdir -p $outdir; fi
# data0=/data/Mimic/version_1_3/processed/mimic_trainData.txt
# meta0=/data/Mimic/version_1_3/processed/mimic_meta.txt
# data=/results/$(basename $data0)
# meta=/results/$(basename $meta0)
# ln - s $data0 $data
# ln - s $meta0 $meta
# $mixehr -f $data -m $meta -i 10000 -k $k -n JCVB0

# plot select topic heatmap (Figure 2a and Figure S4)
Rscript main/Fig2_pheTopic.R
Rscript main/Fig2a_topicHeatmap_sel.R

# plot age-correlated topic (Figure 2b and Figure S5)
Rscript main/Fig2b_ageTopics.R

# plot codes associated with schizo or PTSD (Figure 2c and Figure S6-S14)
Rscript main/Fig2c_phewan_comorbid.R

# plot top patients (Figure 3a and c)
Rscript main/Fig3ac_riskPatients.R

# get select word clouds (Figure 3b)
Rscript main/Fig3b_diseaseTopic_wordclouds.R

# predicting EHR code (Figure 4b Table S3)
Rscript main/Fig4b_impmimic_eval_comparison.R

# predicting prospective mortality (Figure 5a)
Rscript main/Fig5a_predictDeath_summary.R

# predicting prospective mortality (Figure 5b)
Rscript main/Fig5b_predictDeath_diseaseTopics.R

# generate Fig. S1
Rscript model_selections/FigS1_loglik_comparison.R

# generate Fig. S2
Rscript model_selections/FigS2_loglik_comparison_sjcvb0.R

# visualize imputed EHR code (Figure S15)
Rscript impmimic/FigS15_visualize_impmimic.R

# predicting lab results (Figure S16)
Rscript impmimic/FigS16_impmimic_eval_labres.R

# predicting EHR code (Figure S17, S18)
Rscript impmimic/FigS17_impmimic_eval.R
