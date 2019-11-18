# Inferring multimodal latent topics of electronic health records (Li et al., 2019)

This archive reproduces the results of a peer review article.

To execute everything, press 'Run'. This executes `run.sh`, a shell script calling R scripts in sequence. Click the 'Environments and Depdencies' (or 'Env') button to see details of the computational environment.

## Included Materials

In the `/code` pane, in addition to this readme, you will find:

* `mixher`, a directory that contains the C++ source code and Makefile for the MixEHR software;
* `main`, a directory that contains all of the R scripts used to reproduce the main result figures (Figure 1 is a conceptual figure):
    * `Fig2a_topicHeatmap_sel.R` produces Figure 2 panel a and Supplementary Figures S4;
    * `Fig2b_ageTopics.R` produces Figure 2 panel b and Supplementary Figures S5;
    * `Fig2c_phewan_comorbid.R` produces Figure 2 panel c and Supplementary Figures S6-14;
    * `Fig3ac_riskPatients.R` produces Figure 3 panel a and c;
    * `Fig3b_riskPatients.R` produces Figure 3 pane b and Supplementary Figures S3 (when using all of the 75 topics not plotted by default to avoid cluttering the output files);
    * `Fig4b_impmimic_eval_comparison.R` produces Figure 4 panel b (panel a is a conceptual figure);
    * `Fig5a_predictDeath_summary.R` produces Figure 5 panel a
    * `Fig5b_predictDeath_diseaseTopics.R` produces Figure 5 panel b
* `impmimic`, a directory that contains all of the scripts used to perform the full imputation analysis, which eventually lead to Figure 4
    * `createCrossValidSets.R` create 5-fold training and validation sets
    * `select_impute_target_code.R` select relevant EHR target code from each data category for the 5-fold CV
    * `cvimpute.sh` train mixehr on each fold and then infer patient topic and predict missing EHR code on the validation fold (repeat 5 times for the 5-fold CV)
    * `FigS15_visualize_impmimic.R` visualize the predicted EHR code versus the observed EHR code and generates supplementary Figure S15
    * `FigS16_impmimic_eval_labres.R` evalaute missing lab result predictions and generates supplementary Figure S16
    * `FigS17_impmimic_eval.R` evaluate predictions preformance separated by data types and generates Supplementary Figure S17 and S18
    * datatype_flattned, a directory contains the scripts to train the same model but on flattned data for baseline comparison
* `mortality`, a directory that contains all of the scripts to reproduce the mortality prediction results
    * `infer_pat_mix.sh`, trains MixEHR on the training data and then infer patient mixture of test patients using their first admission records
    * `predictDeath_LDA.R`, train LDA model on the flattend EHR data and predict mortality from using their patient topic mixture by LASSO
    * `predictDeath_MXR.R`, predict mortality of test patient using their patient topic mixture by LASSO
* `model_selections`
    * FigS1_loglik_comparison.R, produces Suppl. Figure S1
    * FigS2_loglik_comparison_sjcvb0.R, produces Suppl. Figure S2
* `run.sh`, a [shell script](https://help.codeocean.com/user-manual/whats-up-with-sh-files-on-code-ocean) that the capsule will run by default;

In the `/data` pane, you will find:

* Mimic, a directory contains MIMIC-III data
    * version_1_3/processed
        * `ehrFeatId.RData`, mapping between numerical id of the EHR code and the text description of the code, this is needed to plot the topic definitions learned from MixEHR (e.g., Fig 2a)
        * `mimic_meta.txt`, meta info of the EHR code recognized by the MixEHR software see /code/mixehr/README.txt for details
        * `mimic_trainData.txt`, coded training data recognized by the MixEHR software see /code/mixehr/README.txt for details, these are the records for patients with only one admission
        * `mimic_testData_firstAdm.txt`, EHR code of the first admission for patients with multple admissions
        * `mimic_testData_lastAdm_deathLabel.txt`, death label in the last admission records of the same patients in mimic_testData_firstAdm.txt. This is used to evaluate the model performance
        * `PATIENTS.RData`, detailed information about the patients with single and multiple admissions as well as their ages. The ages are used to plot the age-correlated topics in Fig. 2b
        * `metadata`, a dirctory contains information about the ICD9 code such as group level information (see Fig. S18)
    * ADMISSIONS.csv.gz, detailed admission information about the patient and mortalitiy (i.e., DEATHTIME and HOSPITAL_EXPIRE_FLAG). This file is used in evaluating the mortality predictions (Fig. 5a)
* results, full trained model results
    * mixmimic, directory contains the trained model parameters that are used to generate the topic definitions for downstream analyses
    * impmimic, directory contains cross-validation results for evaluating EHR code predictions in Fig. 4
    * mortality, directory contains the mortality predictions using 50-topic and 75-topic MixEHR and LDA

### Whom to contact:
Any and all questions can be directed to Yue Li at yueli@cs.mcgill.ca.