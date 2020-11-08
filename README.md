---
title: "UniD: Unified Diagnostic Platform of gliomas"
author: "Jie Yang"
---


<a id="top"></a>
[Introduction](#Introduction)

[Installation](#source)

[Data Processing and QC](#QC)

1. [Data loading](#QC1)

2. [Checking Internal Controls](#CheckInternal)

3. [Checking Probe Quality](#CheckQC)

4. [Data normalization (BMIQ)](#BMIQ)

5. [Probe filtering](#filter)

[Data prediction](#prediction)

[Run with Examples](#example)

[System Requirements](#system)



<a id="Introduction"></a>

## Introduction

`UniD` is short for Unified Diagnostic Platform of gliomas, which is designed to predict other genomic information using Illumina Infinium DNA methylation microarray data. This package includes two section of function: data processing QC and genomic information prediction. 

In the data processing and quality control section:

- `UniD.load()` can load the raw data; 
- `UniD.intctl()` can check the quality of internal control probes for each sample; 
- `UniD.dataqc` can evaluate the probe quality in terms of the probe detection P-value and beadcount; 
- `UniD.probefilter()` can filter out probes belong to certain categories and high proportion of missing values and samples with high proportion of missing values; 
- `UniD.BMIQ()` can normalize the data generated from Type I and Type II probes.

In the genomic information prediction section: DNA methylation data generated from both `HumanMethylation 450k` and `HumanMethylation EPIC` platform can be used for prediction. With DNA methylation data, `UniD.pred()` can predict:

- Chromosome 1p/19q co-deletion
- _IDH_ mutation status
- _ATRX_ mutation status
- _TERT_ promoter mutation status
- TCGA gene expression subtype (Classical, Mesenchymal, and Proneural)
- _MGMT_ expression level of tumor

[Back to Top](#top)

<a id="source"></a>
## installation from source file
Then install the package through source file: https://drive.google.com/file/d/13NiEX8SA8LxZ5ayawDpzlAL-HiycxHoD/view?usp=sharing
```
install.packages("/path/to/file/UniD_0.0.1.tar.gz", repos = NULL, type="source")
```
Some prerequisite packages need to be installed from Bioconductor. You can install them such as:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("wateRmelon") ## installing wateRmelon
```
**NOTE**: Java development kit is requried for some prerequisite packages.

<a id="QC"></a>

## Data processing and QC

<a id="QC1"></a>

__1. Data loading__

Raw `.idat` files should be stored in one folder with a `.csv` file which includes the information of sample, their corresponding sentrix ID and sentrix position. This sample sheet should be saved in specific format with `,` as the delimiter. For example (as shown in textEdit):

```
[Header],,,,,,
Investigator Name, Jie Yang,,,,,
Project Name, Example,,,,,
Experiment Name, test,,,,,
Date,01-01-1990,,,,,
[Data],,,,,,
Sample_Name,Sample_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position
sample1,,,,,1001001001,R01C01
sample2,,,,,1001001001,R02C01
sample3,,,,,1001001001,R03C01
sample4,,,,,1001001001,R04C01
```
A list with all necessary information will be created using `UniD.load()`. For example:

```
loading <- UniD.load(dataDir = "~/Desktop/IDAT", outDir = "~/Desktop/output/",
arrayType = "450k", write = T)
```
**note**: the input dataDir should not follow with "/", for example:  
Right version: ~/Desktop/IDAT  
Wrong version: ~/Desktop/IDAT/  

[Back to Top](#top)

<a id="CheckInternal"></a>

__2. Checking Internal Controls__

Different categories of internal control were hold within each beadchip. Based on the criteria provided in Beadchip Control Reporter (BACR), each following categories will be checked and some of them have more than one probes:

- restoration (distinguished by sample Type, FFPE samples needed)
- staining
- extension
- hybridization
- target removal
- bisulfite conversion
- specificity
- nonpolymorphic

However, it is not recommended to excludes sample which failed any categories. However, they could obviously be used to check the experiment quality especially when something goes wrong. A table with all the sample QC for those categories will be returned with function `UniD.intctl()`, as shown in example:

```
samQC <- UniD.intctl(loading = loading, dataDir = "~/Desktop/IDAT/",
outDir = "~/Desktop/output/", arrayType = "450k", sampleType = "FFPE",
write = T)
```

[Back to Top](#top)

<a id="CheckQC"></a>

__3. Checking Probe Quality__

Probes will be checked for detection P-value `detP` and beadcount `bc` in a sample-wise manner which means each sample are independent when being evaluated. 

The detection P-value is calculated by comparing each probes to the general backgroun probes. Probes failed when they have a `detP` > 0.05 which means they are not significant detectable comparing to the background signal. The probe will be set as missing value. The beadcount is the number of bead for for probes. If the `bc` < 3, that specific probe is failed and will be set as missing value for following analysis. 

`Beta.raw` will be generated with the function `UniD.dataqc` with columns as samples and rows as probes. Missing values were included and data was represented using Beta value.

```
Beta.raw <- UniD.dataqc(loading = loading, outDir = "~/Desktop/output/",
detP.cut = 0.05, bc.cut = 3, arrayType = "450k", write = T)
```

[Back to Top](#top)

<a id="BMIQ"></a>

__4. Data normalization (BMIQ)__

Illumina `HumanMethylation 450k` and `HumanMethylation EPIC` platform have both Type I and Type II probes on the microarray while `HumanMethylation 27k`. Due to the chemical differences between the two types of probes exist in the beadchip, Type I and Type II probes have different value distribution. The difference should be OK if all samples were run on the same platform. However, if we want to compare the different samples between `HM450k` and `HM27k`, data should be normalized before any comparison (especially because many samples from TCGA were run on both `HM27k` and `HM450k`).

Many normalization method is available right now, including Beta-mixture quantile normalization `(BMIQ)`, subset-quantile with array normalizaton `SWAN`, peak-based normalization `PBC`, functional normalization, and so on. Because this package is focus on the genomic information prediction, therefore, we only provide one normalization method `BMIQ` here, the gene expression subtype prediction need the specific `BMIQ` normalized data as input. 

This normalization step is recommended to be processed before the probe filtering (`UniD.probefilter`). Because the normalization depends on the beta value distribution of all Type I and Type II probes exist. If the probe filtering processed first, a lot of probes may be filtered out and the normalization process will be easily biased.

```
Beta.BMIQ <- UniD.BMIQ(Beta.raw, outDir = "~/Desktop/output/",
detP.cut = 0.05, bc.cut = 3, arrayType = "450k", write = T)
```

[Back to Top](#top)

<a id="filter"></a>

__5. Probe filtering__

After we generated Beta.raw and Beta.BMIQ, we could use the `UniD.probefilter` to do the probe filtering. Probes can be optionally filtered out if they belonging to the following categories. Each filter is controled by one argument in the function of `UniD.probefilter()`.

- `filterXY`: probes located on Chromosome X or Y (`HM450` and `HM EPIC` different version)
- `filterSNPHit`: probe quality may get affected by SNP nearby (`Zhou W, Nucleic Acids Research, 2017`, `HM450` and `HMEPIC` different version)
- `filterMultiHit`: probes can be mapped to multiple locations (`Nordlund J, Genome Biology, 2013`, `HM450` version)
- `filterNonCG`: probes not targeted to the CpG methylation sites
- `filterNonEpic`: probes available on `HM450` but not available on `HMEPIC`. This is useful when using data from `HM450` to build models and make sure it is suitable for data from `HMEPIC`.

- `filterSample` and `filterSample.cut`: Whether filter out samples with large proportion of missing values (missing values may caused by non-significant detection P-value or low beadcount). The `filterSample.cut` is controling the threshold to filter sample
- `filterProbe` and `filterProbe.cut`: whether filter out probes with large proportion of missing values among samples. This should be use with careful if sample size is small. `filterProbe.cut` is controling the threshold for probe filtering.

```
Beta clean <- UniD.probefilter(Beta.raw, outDir = "~/Desktop/output/",
filterXY = T, filterSNPHit = T, filterMultiHit = T, filterNonCG = F,
filterNonEpic = T, arrayType = "450k", filterSample = T, filterSample.cut
= 0.1, filterProbe = F, write = T)
```

After `UniD.probefilter`, the `Beta.clean` is ready for use as data frame.

[Back to Top](#top)

<a id="prediction"></a>

## Data prediction

All prediction models were compiled within one function `UniD.pred()`. Each argument control the prediction for each genomic information as shown below:

- `Pred.1p19q`: whether predict chromosome 1p/19q co-deletion
- `Pred.IDH`: whether predict _IDH_ mutation status
- `Pred.ATRX`: whether predict _ATRX_ mutation status
- `Pred.TERTp`: whether predict _TERT_ promoter mutation status
- `Pred.ExpressSubtype`: whether predict TCGA gene expression subtype
- `Pred.MGMTExpress`: whether predict _MGMT_ expression level by tumor

However, what need to be emphasize here is we have two different input data sets. For `Pred.ExpressSubtype`, we need to use the `Beta.BMIQ` which is the beta value normalized by `BMIQ` method. For all other predicton model, `Beta.clean` or `Beta.raw` is OK.                       

[Back to Top](#top)

<a id="example"></a>

## Run with examples

The example data set generated from EPIC platform is downloaded from the Illumina website, the demo data set <https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html>. If processing the 450k data, the general data processing workfow is the same as EPIC data set, but need to change the `arrayType` to `450k` and may set the `filterNonEpic` as `TRUE`.

__0. One step function__

For usage convience, all the functions can be called one by one with one function. In this function, all arguments will use the default values. Use with careful.

```
library(UniD)
library_path = .libPaths()
dataDir = paste0(library_path, "/UniD/extdata")

pred <- UniD(dataDir = dataDir, outDir = "./", arrayType = "EPIC", write = T, sampleType = "other")
```

__1. Data loading__

```

loading <- UniD.load(dataDir = dataDir, outDir = "./", arrayType = "EPIC", write = T)

```

    #R log
    ====Data Loading Start====
    Package version of UniD is: 0.0.1
    Loading data from "C:/Program Files/R/R-3.6.2/library/UniD/extdata"
    [read.metharray.sheet] Found the following CSV files:
    [1] "C:/Program Files/R/R-3.6.2/library/UniD/extdata/sample_sheet_EPIC.csv"
    loading(rgSet, Mset, detP) saved as: .//loading.RData
    ====Data Loading Finsihed====
    
__2. Checking Internal Controls__


```
samQC <- UniD.intctl(loading = loading, dataDir = dataDir,
                     outDir = "./", arrayType = "EPIC", sampleType = "other",
                     write = T)
```

    #R log
    ====Internal Control Checking Start====
    Failed sample in Internal Control:

    Restoration_green : 
    staining_green : 
    staining_red : 
    extension_green : 
    extension_red : 
    hybridization_green_HM : 
    hybridization_green_ML : 
    targetRemoval_green1 : 
    targetRemoval_green2 : 
    bisulfite_conversion_I_green1 : 
    bisulfite_conversion_I_green2 : 
    bisulfite_conversion_I_red1 : 
    bisulfite_conversion_I_red2 : 
    bisulfite_conversion_II_1 : 
    bisulfite_conversion_II_2 : 
    specificity_I_green : 
    specificity_I_red : 
    specificity_II_1 : 
    specificity_II_2 : 
    nonpolymorphic_green : 
    nonpolymorphic_red : 

    Internal Controls checking results saved as: .//samQC_Internal.Control.csv
    ====Internal Control Checking Finished====


__3. Checking Probe Quality__

```
Beta.raw <- UniD.dataqc(loading = loading, outDir = "./",
                        detP.cut = 0.05, bc.cut = 3, arrayType = "EPIC", write = T)
```

    #R log
    ====Sample-wise Probe Quality Assessment Start====
    [1] "The failed fraction per sample (failed detP and bc may overlap): "
                        Fail.Frac.detP Fail.Frac.beadcount Fail.Frac.NA
    200144450018_R04C01    0.001067099         0.001690054  0.002669478
    200144450019_R07C01    0.001085557         0.001621991  0.002609490
    200144450021_R05C01    0.000711784         0.001390113  0.002070749
    The raw Beta value were saved as: .//UniD_Beta.raw.RData
    The fraction of failed probes per sample saved as: .//Failed_probe_Fraction.csv
    ====Sample-wise Probe Quality Assessment Finished====

__4. Data normalization (BMIQ)__

```
Beta.BMIQ <- UniD.BMIQ(Beta.raw, outDir = "./",
                        arrayType = "EPIC", write = T)
```


    #R log
    ====BMIQ Normalization Start====
    [1] "BMIQ on sample: 200144450018_R04C01"
    200144450018_R04C01 has missing value: 2314
    [1] "Fitting EM beta mixture to type1 probes"
    [1] Inf
    [1] 0.005296374
    ... 
    ...
    [1] "Done"
    [1] "Start normalising type 2 probes"
    [1] "Finished for sample 200144450021_R05C01"
    BMIQ normalized Beta value saved as: .//UniD_Beta.BMIQ.RData
    ====BMIQ Normalization Finished====


__5. Probe filtering__

```
Beta.clean <- UniD.probefilter(Beta.raw = Beta.raw, outDir = "./",
                              filterXY = T, filterSNPHit = T, filterMultiHit = T,
                              filterNonCG = F, filterNonEpic = F, arrayType = "EPIC",
                              filterSample = T, filterSample.cut = 0.1, filterProbe = F,
                              write = T)
```

    #R log
    ====Probe Filter Start====
    Delete 19681 probes on ChrX/Y for EPIC Platform.
    Delete 78720 probes affected by SNP for EPIC Platform. (Zhou W, Nucleic Acids Research, 2017)
    Delete 63 probes mapped to multiple site for EPIC Platform. (Nordlund J, Genome Biology, 2013)
    [1] "The failed fraction per sample (after all probe filters): "
                        Fail.Frac.NA
    200144450018_R04C01  0.002178632
    200144450019_R07C01  0.002108354
    200144450021_R05C01  0.001668463

    Final data matrix is: 3 samples * 768372 probes.
    ====Probe Filter Finish====


__6. Data prediction__

```
Pred <- UniD.pred(inputdata = Beta.raw, inputdata.BMIQ = Beta.BMIQ, inputvalueType = "B",
                  Pred.IDH = T, Pred.1p19q = T, Pred.ATRX = T, Pred.TERTp = T,
                  Pred.ExpressSubtype = T, 
                  outDir = "./", write = T)
```

    #R log
    ====Biomarker Prediction Start====
    #Note: inputdata and inputdata.BMIQ must with columns as samples
              and rows as probes
    #Note: inputdata.BMIQ must be Beta value format

    ##Predict Chromosome 1p19q co-deletion Start##
    Change B value to M value
    No missing value for predictor.
    #Predict Chromosome 1p19q co-deletion Finish##

    ##Predict IDH mutation status Start##
    Change B value to M value
    KNN used to impute 1 missing value
    Imputation Done (Check missing value fraction for each sample)
    ##Predict IDH mutation status Finish##

    ##Predict ATRX mutation status Start##
    Change B value to M value
    KNN used to impute 4 missing value
    Imputation Done (Check missing value fraction for each sample)
    ##Predict ATRX mutation status Finish##

    ##Predict TERT promoter mutation status Start##
    Change B value to M value
    KNN used to impute 11 missing value
    Imputation Done (Check missing value fraction for each sample)
    ##Predict TERT promoter mutation status Finish##

    ##Predict TCGA Gene Expression subtype Start##
    Change B value to M value
    No missing value for predictor.
    The following `from` values were not present in `x`: 1, 3
    ##Predict TCGA Gene Expression subtype Finish##

    Predicted result was saved as: .//UniD_Biomarker.Pred.csv
    ====Biomarker Prediction Finish====


<a id="system"></a>
## System requirements

### OS requirements
Our UniD package has been tested on the following environment:

- macOS: Mojave 10.14.6
- Windows 10. R version 3.6.2

### Dependent R packages
Please see the DESCRIPTION.

[Back to Top](#top)

