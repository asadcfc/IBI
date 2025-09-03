#!/usr/bin/env python
# coding: utf-8

# In[1]:


### IBI implementation in python: Developed by Jinling Liu 12/26/2021
### Can be applied to multiple traits or single trait
### 02-09-22 updated the function of lgM_cal_1 so it can either calculate using topGD or sGD, with a flag of "use_topGD" before that function


# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn2_circles
import json
import os
import shutil
import scipy.stats as stats
from scipy.stats import fisher_exact
import scipy
from datetime import datetime
import math
import time
from joblib import Parallel, delayed


# In[3]:


def read_traitsF (traitsF):  ### read the .csv file with traits (subjects_row x traits_column)
    traits = pd.read_csv(traitsF,index_col = 0) 
    subIDs = list(traits.index)
    traitIDs = traits.columns
    traits = np.array(traits,dtype=np.int16)
    
    return (subIDs, traitIDs, traits) 


# In[4]:


def read_variantsF1 (variantsF):## read the large genomic file (row_SNPs x column_subjects) using pandas
    df = pd.read_csv(variantsF,index_col = 0) 
    varIDs = list(df.index)
    subIDs = list(int(x) for x in df.columns)
    variants = np.array(df,dtype=np.int8) 
    A0 = np.ones(len(subIDs),dtype=np.int8)
    variants = np.row_stack((A0,variants))
    varIDs.insert(0,'A0')
    df = pd.DataFrame(variants,index=varIDs,columns=subIDs,dtype=np.int8)
    
    return (subIDs, varIDs, variants, df) 


# In[5]:


def DriverSearch (traits,variants): 
    ### Calcuate and return the lgM for all the drivers or any driver for any given population for multiple traits
    bpMask0 = traits==0
    bpMask0 = bpMask0.astype(np.int16)
    d0 = np.sum(bpMask0)
    bpMask1 = traits==1
    bpMask1 = bpMask1.astype(np.int16)
    d1 = np.sum(bpMask1)
    snpMask0 = variants==0
    snpMask0 = snpMask0.astype(np.int16)
    snpMask1 = variants==1
    snpMask1 = snpMask1.astype(np.int16)
    V0D0 = snpMask0@bpMask0
    V1D0 = snpMask1@bpMask0 
    V0D1 = snpMask0@bpMask1
    V1D1 = snpMask1@bpMask1
    # when j=0 (V=0)
    lgM = scipy.special.loggamma(2.0) - scipy.special.loggamma(2.0+V0D1+V0D0)
    lgM += scipy.special.loggamma(1.0+V0D0) - scipy.special.loggamma(1.0)
    lgM += scipy.special.loggamma(1.0+V0D1) - scipy.special.loggamma(1.0)
    # when j=1 (V=1)
    lgM += scipy.special.loggamma(2.0) - scipy.special.loggamma(2.0+V1D1+V1D0)
    lgM += scipy.special.loggamma(1.0+V1D0) - scipy.special.loggamma(1.0)
    lgM += scipy.special.loggamma(1.0+V1D1) - scipy.special.loggamma(1.0)
    if variants.ndim == 1:
        lgM = lgM.reshape(1,lgM.shape[0])
        
    return(lgM)


# In[6]:


def GDsearch_all(traits,variants): 
    ## Get all the stats for all the variants in any given population for multiple traits; particulary used for the entire population
    bpMask0 = traits==0
    bpMask0 = bpMask0.astype(np.int16)
    d0 = np.sum(bpMask0)
    bpMask1 = traits==1
    bpMask1 = bpMask1.astype(np.int16)
    d1 = np.sum(bpMask1)
    snpMask0 = variants==0
    snpMask0 = snpMask0.astype(np.int16)
    snpMask1 = variants==1
    snpMask1 = snpMask1.astype(np.int16)
    V0D0 = snpMask0@bpMask0
    V1D0 = snpMask1@bpMask0
    V0D1 = snpMask0@bpMask1
    V1D1 = snpMask1@bpMask1
    cp_D1V1 = (1 + V1D1)/(2 + V1D1 + V1D0)*1.0                 
    cp_D1V0 = (1 + V0D1)/(2 + V0D1 + V0D0)*1.0                
    RR = cp_D1V1/cp_D1V0 
    # when j=0 (V=0)
    lgM = scipy.special.loggamma(2.0) - scipy.special.loggamma(2.0+V0D1+V0D0)
    lgM += scipy.special.loggamma(1.0+V0D0) - scipy.special.loggamma(1.0)
    lgM += scipy.special.loggamma(1.0+V0D1) - scipy.special.loggamma(1.0)
    # when j=1 (V=1)
    lgM += scipy.special.loggamma(2.0) - scipy.special.loggamma(2.0+V1D1+V1D0)
    lgM += scipy.special.loggamma(1.0+V1D0) - scipy.special.loggamma(1.0)
    lgM += scipy.special.loggamma(1.0+V1D1) - scipy.special.loggamma(1.0)
    if variants.ndim == 1:
        lgM = lgM.reshape(1,lgM.shape[0])
    max_value = np.max(lgM,axis=0) 
    max_index = np.argmax(lgM,axis=0)
    
    return (RR,lgM,max_value,max_index)


# In[7]:


def lgMcal_1(varID): ## use DriverSearch for V0 and SD_lgM_V1 for V1 ## designed for using one topGD
    i = varIDs.index(varID) 
    index1 = variants[i,:]==1
    index0 = variants[i,:]==0
    V1 = variants[:,index1]
    if use_oneTopGD:
        V0 = variants[topGD_index][:,index0] 
    else:
        V0 = variants[:,index0] 
    BP_V1 = traits[index1] 
    BP_V0 = traits[index0]
    lgMv1_SD = DriverSearch(BP_V1,variants[i,index1])[0]
    lgMv0 = DriverSearch(BP_V0,V0) 
    lgMv0_topGD = [] 
    r = []           
    if use_oneTopGD:
        for m in range(0,len(traitIDs)):
            lgMv0_topGD.append(lgMv0[m,m])
        for j in topGD_index: 
            r1 = stats.spearmanr(variants[i,:],variants[j,:])[0]
            r.append(r1)
        lgMv0_sGD = np.zeros(len(traitIDs))
        sGD = np.zeros(len(traitIDs))
    else:
        lgMv0_sGD = np.max(lgMv0,axis=0) 
        sGD_index = np.argmax(lgMv0,axis=0)          
        sGD = [] 
        for item in sGD_index:
            sGD.append(varIDs[item])
        sGD = np.array(sGD) 
        k=0  
        for j in topGD_index: 
            lgMv0_topGD.append(lgMv0[j,k]) 
            r1 = stats.spearmanr(variants[i,:],variants[j,:])[0]
            r.append(r1)
            k = k + 1 
    lgMv0_topGD = np.array(lgMv0_topGD) 
    r = np.array(r) 
    lgM_v1v0 = lgMv1_SD + lgMv0_sGD
    
    return(lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID)


# In[ ]:


## Main file: use the above functions to get all statistics

subIDs, varIDs, variants, df_variants = read_variantsF1 ('file_path/Example_Variants_file_name.csv')
subIDs_BP, traitIDs, traits = read_traitsF('file_path/Example_Traits_file_name.csv')

rr, glgm, glgm_topGD, topGD_index = GDsearch_all(traits,variants)

use_oneTopGD = True
lgMv1_SD, lgMv0_sGD, lgMv0_topGD, lgM_v1v0, sGD, r, i, varID = lgMcal_1('varID')

element_run = Parallel(n_jobs=-1)(delayed(lgMcal_1)(var) for var in varIDs[0:100]) 

