# ROMA
Representation of Module Activity.

## Introduction ##

ROMA is a software package written in Java for the quantification and representation of biological module activity using expression data.

In many biological studies, there is a need to quantify the activity of a set of genes in individual samples. For instance, we often would like to quantify the effect of a given regulator, e.g. an active transcription factor, on a pool of target genes. 

ROMA is using the first principal component of a PCA analysis to summarize the coexpression of a group of genes in the gene set. ROMA is also providing additional functionalities: (i) calculation of the individual gene contribution to the module activity level and determination of the genes that are contributing the most to the first pricipal component, (ii) Alternative techniques to compute the first principal component, i.e. weighted and centered methods, (iii) estimation of the statistical significance of the proportion of variance explained by the first principal component, as well as the spectral gap between the variance explained by the first and second component (representing the homogeneity for the gene set).

## Compilation and installation ##

First clone the GitHub repository on your computer with the command:

```
git clone https://github.com/sysbio-curie/Roma/ .
```

The directory contains a configuration file for building with the [ant](https://en.wikipedia.org/wiki/Apache_Ant) tool. You can compile the files 
