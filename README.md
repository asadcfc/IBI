# IBI
Individualized Bayesian Inference (IBI) algorithm is used for estimating the genomic variants that influence complex traits such as hypertension at the level of an individual (e.g., a patient). It is based on Python libraries.

# Example
subIDs, varIDs, variants, df_variants = read_variantsF1 ('file_path/Example_Variants_file_name.csv')
subIDs_BP, traitIDs, traits = read_traitsF('file_path/Example_Traits_file_name.csv')

rr, glgm, glgm_topGD, topGD_index = GDsearch_all(traits,variants) 

lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID = lgMcal('varID')

use_oneTopGD = True
lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID = lgMcal_1('varID')
