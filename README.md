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

The directory contains a configuration file for building with the [ant](https://en.wikipedia.org/wiki/Apache_Ant) tool. You can compile the files and build the jar file at the same time with:

```
ant jar
```

If everything went well, you should see a roma_v1.0.jar file in the directory.

You can now run the program with:

```
java -jar roma_v1.0.jar
```

This should give you a list of the options for ROMA:

```
:: OPTIONS
REQUIRED:
:: -dataFile : data in tab-delimited format (required)
:: -moduleFile : name of the gmt module file  (required)
Optional:
:: -outputFolder : folder name for keeping the resulting files (by default it will be the folder of the data file)
:: -sampleFile : description of samples in tab-delimited txt format
:: -saveDecomposedFiles : save dat files for each module in the output folder
:: -outlierThreshold : threshold for determining outliers
:: -typeOfModuleFile : 0 - for standard GMT, 1 - for GMT with weights (default)
:: -typeOfPCAUsage : 0 - for standard PCA, 1 - for PCA with fixed center (default)
:: -robustPCA : 0 - all points are used, 1 - leave one out -based removal of outliers
:: -robustPCASampling : 0 - all points are used in random sampling, 1 - leave one out -based removal of outliers in random sampling (slow down calculations)
:: -centerData : 0 - do not center, 1 - center each line
:: -mostContributingGenesZthreshold : (default 1) threshold (z-value) used to print out the names of the genes most contributing to the component
:: -diffSpotGenesZthreshold : (default 1) threshold (t-test) used to print out the names of the differentially expressed genes
:: -fieldForAveraging : (optional) number of the field column used for computing the average module activities (ex: 5)
:: -fieldValueForAveraging : (optional) value of the field used for computing the reference value of module activity (ex: "normal")
:: -fieldForDiffAnalysis : (optional) number of the field column used for differential analysis (ex: 5)
:: -fieldValuesForDiffAnalysis : (optional) values of the field column used for differential analysis separated by % symbol (ex: "invasive%noninvasive")
:: -numberOfPermutations : (optional) number of samples for empirical p-values estimation
:: -numberOfGeneSetSizesToSample : (optional) number of random gene set sizes to test for empirical p-values estimation
:: -minimalNumberOfGenesInModule: minimal size of the gene set
:: -minimalNumberOfGenesInModuleFound: minimal number of genes in a gene set found in dataset

```

If you want to run the program fro another directory, you can create a link to the roma jar file, and run it with the same command in the distant directory. 

Another way would be to set the CLASSPATH to include the roma jar file and necessary libraries:

```
export CLASSPATH="/Users/johndoe/git_roma/roma_v1.0.jar:/Users/johndoe/git_roma/lib/VDAOEngine.jar"
```

Then you should be able to run ROMA anywhere with:

```
java fr.curie.ROMA.ModuleActivityAnalysis
```

