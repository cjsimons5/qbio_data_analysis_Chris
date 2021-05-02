#!/usr/bin/env python
# coding: utf-8

#Spearman Rho Correlation Plotting Function for Specific Gene
def spear_rho_plot(rna, protein, genename, figsz=10, pathout):

    rho = stats.spearmanr(rna, protein)
    #Calculates spearman rho to use given rna and proteomic data for a specific gene
    m, b = np.polyfit(protein, rna, 1)
    #Calculates slope of data points in the correlation
    
    plt.figure(figsize=(figsz, figsz))
    
    #Plots the data points with a green trendline going through it
    plt.scatter(protein, rna, c="black")
    plt.plot(protein, m*protein + b, 'g')
    
    #Graph labels
    title = "rho: {0} for {1}".format(rho_old, genename) #Title of graph using string formatting
    plt.title(title)
    plt.xlabel("Proteomic Data")
    plt.ylabel("RNA Data")
    
    #Saves figure to the path given in function if running as python script
    plt.savefig(pathout, bbox_inches="tight" )

