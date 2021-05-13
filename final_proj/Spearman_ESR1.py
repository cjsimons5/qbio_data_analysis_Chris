#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('jupyter nbconvert --to script Spearman_ESR1.ipynb')


# In[1]:


import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats 


# In[2]:


cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

#The dataframes are MultIndex pandas dataframes. 
#However, to teach the basics of pandas, we will remove the "multi" part of the dataframe.
protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()


# In[3]:

#Gets age data
clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12


# In[4]:


assert list(rna_data.index) == list(protein_data.index)


# In[5]:

#shrinking down data to only include gene ESR1
rna_esr1 = rna_data.loc[: , "ESR1"]
protein_esr1 = protein_data.loc[: , "ESR1"]

rho, spear_pvalue = stats.spearmanr( rna_esr1, protein_esr1 )

rho_check, spear_pvalue_check = stats.spearmanr( protein_esr1, rna_esr1 )

assert rho == rho_check


# In[8]:


#Boolean masking, using age data
young_mask = clinical_data["Age_in_years"] < 50.0
old_mask = clinical_data["Age_in_years"] >= 50.0

#Splitting into young using the mask (at row ESR1, it is only where the mask is true)
rna_esr1_young = rna_data["ESR1"][ young_mask ]
protein_esr1_young = protein_data["ESR1"][ young_mask ]

#Splitting into old using the mask (at row ESR1, it is only where the mask is true)
rna_esr1_old = rna_data["ESR1"][ old_mask ]
protein_esr1_old = protein_data["ESR1"][ old_mask ]


# In[14]:


def spear_rho_plot(rna, protein, genename, pathout, figsz=10):
    #input takes in rna expression and protein expression data, a gene name
    #a pathway to save the plot to, and a figure size with default dimension 10

    #Calculation of the spearman rho
    rho, spear_pval = stats.spearmanr(rna, protein)
    #prepping linear fit line
    m, b = np.polyfit(protein, rna, 1)
    plt.figure(figsize=(figsz, figsz))
    
    #plotting the data points, and plotting the line of best fit
    plt.scatter(protein, rna, c="black")
    plt.plot(protein, m*protein + b, 'g')
    
    title = "rho: {0} for {1}".format(rho, genename) #Puts rho correlation and genename in title of figure
    plt.title(title)

    #informative labels
    plt.xlabel("Proteomic Data")
    plt.ylabel("RNA Data")
    
    #Saves to the path put in
    plt.savefig(pathout, bbox_inches="tight" )


# In[17]:

out="/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphALL.png"
spear_rho_plot(rna_esr1, protein_esr1, "ESR1",  out)

out="/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphOLD.png"
spear_rho_plot(rna_esr1_old, protein_esr1_old, "ESR1",  out)

out="/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphYOUNG.png"
spear_rho_plot(rna_esr1_young, protein_esr1_young, "ESR1",  out)
