#!/usr/bin/env python
# coding: utf-8

# In[3]:


import cptac
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats 


# In[4]:


cptac.download(dataset="Brca")
br = cptac.Brca()

protein_data = br.get_proteomics()

#The dataframes are MultIndex pandas dataframes. 
#However, to teach the basics of pandas, we will remove the "multi" part of the dataframe.
protein_data = protein_data.droplevel(1, axis=1)

rna_data = br.get_transcriptomics()
clinical_data = br.get_clinical()


# In[5]:


clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12


# In[6]:


assert list(rna_data.index) == list(protein_data.index)


# In[21]:


rna_esr1 = rna_data.loc[: , "ESR1"]
protein_esr1 = protein_data.loc[: , "ESR1"]

rho, spear_pvalue = stats.spearmanr( rna_esr1, protein_esr1 )

rho_check, spear_pvalue_check = stats.spearmanr( protein_esr1, rna_esr1 )

assert rho == rho_check


# In[27]:


plt.figure( figsize=(10,10) )

#Replace x and y with appropriate variables
plt.scatter( protein_esr1, rna_esr1 )

title = "rho: {} for ESR1 (all ages)".format(rho) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

#Fill in informative x and y labels
plt.xlabel("Proteomic Data")
plt.ylabel("RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphALL.png", bbox_inches="tight" )


# In[17]:


#What column of clinical_data is referring to age?
young_mask = clinical_data["Age_in_years"] < 50.0
old_mask = clinical_data["Age_in_years"] >= 50.0

#Check for understanding: Why do the below lines work?
rna_esr1_young = rna_data["ESR1"][ young_mask ]
protein_esr1_young = protein_data["ESR1"][ young_mask ]

#We want all patients of the ESR1 column
rna_esr1_old = rna_data["ESR1"][ old_mask ]
protein_esr1_old = protein_data["ESR1"][ old_mask ]


# In[23]:


#YOUNG PLOT
rho_young, spear_pvalue_young = stats.spearmanr( rna_esr1_young, protein_esr1_young )

plt.figure( figsize=(10,10) )

#Replace x and y with appropriate variables
plt.scatter( protein_esr1_young, rna_esr1_young )

title = "rho: {} for ESR1 (Patients < 50 years old)".format(rho_young) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

plt.xlabel("Young Proteomic Data")
plt.ylabel("Young RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphYOUNG.png", bbox_inches="tight" )


# In[24]:


#OLD PLOT
rho_old, spear_pvalue_old = stats.spearmanr( rna_esr1_old, protein_esr1_old )

plt.figure( figsize=(10,10) )

#Replace x and y with appropriate variables
plt.scatter( protein_esr1_old, rna_esr1_old )

title = "rho: {} for ESR1 (Patients >= 50 years old)".format(rho_old) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

plt.xlabel("Old Proteomic Data")
plt.ylabel("Old RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphOLD.png", bbox_inches="tight" )

