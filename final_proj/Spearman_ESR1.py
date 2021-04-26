#!/usr/bin/env python
# coding: utf-8

# In[28]:


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


clinical_data["Age_in_years"] = clinical_data["Age.in.Month"]/12


# In[4]:


assert list(rna_data.index) == list(protein_data.index)


# In[5]:


rna_esr1 = rna_data.loc[: , "ESR1"]
protein_esr1 = protein_data.loc[: , "ESR1"]

rho, spear_pvalue = stats.spearmanr( rna_esr1, protein_esr1 )

rho_check, spear_pvalue_check = stats.spearmanr( protein_esr1, rna_esr1 )

assert rho == rho_check


# In[10]:


plt.figure( figsize=(10,10) )

m, b = np.polyfit(protein_esr1, rna_esr1, 1)

#Replace x and y with appropriate variables
plt.scatter( protein_esr1, rna_esr1, c='black')
plt.plot(protein_esr1, m*protein_esr1 + b, 'g')

title = "rho: {} for ESR1 (all ages)".format(rho) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

#Fill in informative x and y labels
plt.xlabel("Proteomic Data")
plt.ylabel("RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphALL.png", bbox_inches="tight" )


# In[11]:


#What column of clinical_data is referring to age?
young_mask = clinical_data["Age_in_years"] < 50.0
old_mask = clinical_data["Age_in_years"] >= 50.0

#Check for understanding: Why do the below lines work?
rna_esr1_young = rna_data["ESR1"][ young_mask ]
protein_esr1_young = protein_data["ESR1"][ young_mask ]

#We want all patients of the ESR1 column
rna_esr1_old = rna_data["ESR1"][ old_mask ]
protein_esr1_old = protein_data["ESR1"][ old_mask ]


# In[12]:


#YOUNG PLOT
rho_young, spear_pvalue_young = stats.spearmanr( rna_esr1_young, protein_esr1_young )

plt.figure( figsize=(10,10) )
m, b = np.polyfit(protein_esr1_young, rna_esr1_young, 1)

#Replace x and y with appropriate variables
plt.scatter( protein_esr1_young, rna_esr1_young, c='black' )
plt.plot(protein_esr1_young, m*protein_esr1_young + b, 'g')

title = "rho: {} for ESR1 (Patients < 50 years old)".format(rho_young) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

plt.xlabel("Young Proteomic Data")
plt.ylabel("Young RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphYOUNG.png", bbox_inches="tight" )


# In[13]:


#OLD PLOT
rho_old, spear_pvalue_old = stats.spearmanr( rna_esr1_old, protein_esr1_old )

plt.figure( figsize=(10,10) )
m, b = np.polyfit(protein_esr1_old, rna_esr1_old, 1)

#Replace x and y with appropriate variables
plt.scatter( protein_esr1_old, rna_esr1_old, c='black' )
plt.plot(protein_esr1_old, m*protein_esr1_old + b, 'g')

title = "rho: {} for ESR1 (Patients >= 50 years old)".format(rho_old) #This is string formatting. The variable in the () will print in the {}
plt.title(title)

plt.xlabel("Old Proteomic Data")
plt.ylabel("Old RNA Data")

#plt.show() #Comment out when running in script
plt.savefig( "/Users/Christopher/Desktop/Datanalysis/qbio_data_analysis_Chris/final_proj/SpearGraphOLD.png", bbox_inches="tight" )

