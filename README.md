# feature_selection_bathing_water
feature selection comparison for bathing water quality prediction

This repo contributes to Digital Water City Project from Kompetenzzentrum Wasser Berlin gGmbH.

It compares different feature selection (FS) algorithms in their suitability to select good features for bathing water quality prediction with 6 different datasets in germany.

Tested FS algorithms: 
- four forward stepwise selection approaches. Used criteria: AIC and BIC in two variations each. Limited to maximum 5 features and unlimited
- two lasso models, lasso_min and lasso_1se. Lasso_1se comming form the one standard error rule.
