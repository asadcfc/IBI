# IBI
Individualized Bayesian Inference (IBI) algorithm is used for estimating the genomic variants that influence complex traits such as hypertension at the level of an individual (e.g., a patient). It is based on Python libraries.

# Example
Jupyter-notebook file define the functions. We can use the functions to get the statistics. Some examples are given below

subIDs, varIDs, variants, df_variants = read_variantsF1 ('file_path/Example_Variants_file_name.csv')

subIDs_BP, traitIDs, traits = read_traitsF('file_path/Example_Traits_file_name.csv')

rr, glgm, glgm_topGD, topGD_index = GDsearch_all(traits,variants) 

lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID = lgMcal('varID')

use_oneTopGD = True

lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID = lgMcal_1('varID')

element_run = Parallel(n_jobs=-1)(delayed(lgMcal_1)(var) for var in varIDs[0:100]) 

# Figure 2
Figure 2. A Bayesian method for GWAS analysis. "Figure 2.ipynb" jupyter-notebook file  generate the Figure 2a, 2b, and 2c

# Figure 3
Figure 3. Comparison of IBI and GWAS. "Figure 3.ipynb" jupyter-notebook file  generate the Figure 3c, and 3d.
For Figures 3a and 3b, follow the steps below.:
  Step 1: Use OriginPro software 
  Step 2: Use its Manhattan plot app 
  Step 3: Generate plot for GWAS, and IBI

# Figure 4
Figure 4. HTN patient coverage. "Figure 4.ipynb" jupyter-notebook file  generate the Figure 4a, 4b, and 4c.

# Figure 5
Figure 5. ROC curves for predicting hypertension from top SNPs. "Figure 5.ipynb" jupyter-notebook file  generate the Figure 5a, 5b, 5c, 5d, 5e and 5f.
