# Log-Linear model for biases in microbiome studies
---
Ni Zhao (nzhao10@jhu.edu), Glen Satten (GSatten@emory.edu)
---


## Overview
The MicroBias package contains functions that use log-linear model to assess for bias in microbiome sequencing studies. Currently, the model works when the true relative abundances (including which taxa are present in samples) are known, therefore, is only applicable to mock community data. Extensions to situations when true relative abundances are unknown are under development. Please refer to [1] for details about the methodology. Here, we will use the famous Brooks dataset [2] to demonstrate the utility of this package. 

## Depends
MicroBias depends on package limSolve

## Package download and installation 

We can install the package through github. This install_github functino requires devtools pacakage in R.  Then we load the package as well as limSolve, which this package depends on. 

```{r download, eval = F}
library(devtools)
install_github("zhaoni153/MicroBias")
library(limSolve)
library(MicroBias)
```

## MicroBias model description
The MicroBias is a log-linear model for bias inference in microbiome studies. It is an extension of the multiplicative bias generation model from [2]. The model is specified as 

<img src="https://latex.codecogs.com/gif.latex?\ln{\tilde&space;p_{ij}}&space;=&space;\ln{p_{ij}}&space;&plus;&space;X_{i\cdot}&space;\beta_{\cdot&space;j}&space;&plus;&space;\epsilon_{ij},&space;\text{&space;for&space;}&space;\Delta_{ij}&space;>&space;0," title="\ln{\tilde p_{ij}} = \ln{p_{ij}} + X_{i\cdot} \beta_{\cdot j} + \epsilon_{ij}, \text{ for } \Delta_{ij} > 0," />


in which <img src="https://latex.codecogs.com/gif.latex?\tilde&space;p_{ij}" title="\tilde p_{ij}" /> and <img src="https://latex.codecogs.com/gif.latex?p_{ij}" title="p_{ij}" /> are the observed and the true relative abundance level in sample <img src="https://latex.codecogs.com/gif.latex?i" title="i" /> for taxon <img src="https://latex.codecogs.com/gif.latex?j" title="j" /> (<img src="https://latex.codecogs.com/gif.latex?i&space;\in&space;(1,&space;\cdots,&space;n),&space;j&space;\in&space;(1,&space;\cdots,&space;J)" title="i \in (1, \cdots, n), j \in (1, \cdots, J)" />). <img src="https://latex.codecogs.com/gif.latex?\Delta" title="\Delta" /> is the indicator matrix that <img src="https://latex.codecogs.com/gif.latex?\Delta_{ij}&space;=&space;1" title="\Delta_{ij} = 1" /> is the j-th taxon is present in sample <img src="https://latex.codecogs.com/gif.latex?i" title="i" /> by design, and <img src="https://latex.codecogs.com/gif.latex?\Delta_{ij}&space;=&space;0" title="\Delta_{ij} = 0" /> otherwise. For all samples, <img src="https://latex.codecogs.com/gif.latex?\sum_j{\tilde&space;p_{ij}\Delta_{ij}}&space;=&space;\sum_j{p_{ij}\Delta_{ij}}&space;=&space;1" title="\sum_j{\tilde p_{ij}\Delta_{ij}} = \sum_j{p_{ij}\Delta_{ij}} = 1" /> because of the compositionality nature of the data. Here, we assume that <img src="https://latex.codecogs.com/gif.latex?p_{ij}" title="p_{ij}" /> and <img src="https://latex.codecogs.com/gif.latex?\Delta_{ij}" title="\Delta_{ij}" /> are known, as in mock community studies. <img src="https://latex.codecogs.com/gif.latex?X_{i\cdot}" title="X_{i\cdot}" /> represents sample-level covariates, and <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /> are corresponding coefficients. <img src="https://latex.codecogs.com/gif.latex?\epsilon" title="\epsilon" /> is the random errors with mean zero and finite second moments. 

Suppose the design matrix <img src="https://latex.codecogs.com/gif.latex?X" title="X" /> is a <img src="https://latex.codecogs.com/gif.latex?n&space;\times&space;M" title="n \times M" /> matrix (<img src="https://latex.codecogs.com/gif.latex?n" title="n" /> is the sample size,  <img src="https://latex.codecogs.com/gif.latex?M" title="M" /> is the number of covariates), then <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /> matrix is a <img src="https://latex.codecogs.com/gif.latex?M&space;\times&space;J" title="M \times J" /> matrix with <img src="https://latex.codecogs.com/gif.latex?J" title="J" /> being the number of taxa in the system. Appropriate contrasts of the <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /> can be tested for questions about the bias factors. We will use the Brooks data [2] to demonstrate the use of this package. 


## Brooks data set

Brooks et al. [3] conducted a large-scale model community study using seven bacteria commonly found in vaginal samples, i.e., L crispatus,  L iners, G vaginalis, A vaginae, P bivia, S amnii and group B streptococcus (GBS). Each mock sample was designed to consist a mixture of bacterial taxa. Every sample was processed in three experiments using a common experimental workflow, but starting from different  types of samples: even mixtures of cells, of extracted DNA, and of PCR product. There were totally six plates in which the samples were processed in. Plates 1 \& 2 are cell samples, plates 3 \& 4 are DNA samples and plates 5 \& 6 are PCR product. The study was well-balanced that the samples on all plates are comparable with respect to the compositions of bacteria present.

We included the Brooks data in the package, which can be loaded as following: 

```{r loadBrooks}
data(brooks)
names(brooks)
```

As shown, the brooks dataset is organized into a list with four elements. For all of them, each row represents one sample with a total of 240 samples. The rows are matched so that they are in the same order. 

```{r loadBrooks2}
delta = brooks$delta
meta.data = brooks$meta.data
obs.counts = brooks$obs.counts
true.freq = brooks$true.freq
```

### Data clean up
Sometimes a taxon can have non-zero read even when it was not designed to be present. Further, in the obs.counts dataset, there is a column called "other", which we will remove as well.

```{r cleanData}
otu.tab = obs.counts[,-which(colnames(obs.counts) == "Other")]
colnames(otu.tab)=c('Avaginae', 'Gvaginalis', 'Lcrispatus', 'Liners', 'Pbivia', 'Samnii', 'GroupBStrep')
# Remove taxa that are not present by design
otu.tab = otu.tab * delta

## Calculate relative abundancies
otu.tab=otu.tab/rowSums(otu.tab)

```

Further, we should remove the samples that have only a single taxon present as they contribute no information to the algorithm. 

```{r rmSingle}
multi.otu=which( rowSums(delta)>1 )
otu.tab=otu.tab[multi.otu,]
true.freq=true.freq[multi.otu,]
meta.data=meta.data[multi.otu,]
delta=delta[multi.otu,]
k.s=rowSums(delta)
```

## Least square estimation 
Here, we fit the full model with both plate effect and the taxon-taxon interaction effect. The design matrix has a total of 12 columns, with the first seven representing the <img src="https://latex.codecogs.com/gif.latex?\Delta" title="\Delta" />, and the last 5 columns representing the plate indicator for plates 2 to 6. Indicator for the first plate is not included. The <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /> matrix is a <img src="https://latex.codecogs.com/gif.latex?12&space;\times&space;7" title="12 \times 7" /> matrix with the off-diagonal elements in the first seven rows representing the taxon-taxon interaction, the diagonal elements in the first seven rows and the 8-th to 12-th rows representing the main plate effect.


```{r estimate}
x =cbind(delta , model.matrix(~factor(meta.data$Plate)+0)[,2:6] )
mod = LLBias.fit(otu.tab = otu.tab,x = x, offset = true.freq , delta = delta)
```

## Hypothesis testing

Contrasts of <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /> may be tested if it satisfies the testability requirement, as shown in [1]. Constructing the contrast matrix is usually challenging, as the number of parameters is large in our model. We provide two ways to construct the contrast matrix. 

As an example, we test whether A vaginae interact with the bias of any of the other OTUs. This corresponds to the H7(1) in the simulations of the paper [1]. Recall that our <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /> matrix is as follows: 

<img src="https://latex.codecogs.com/gif.latex?\begin{pmatrix}&space;\beta_1&space;&&space;\beta_{1,2}&space;&&space;\cdots&space;&&space;\beta_{1,7}\\&space;\beta_{2,1}&space;&&space;\beta_2&space;&&space;\cdots&space;&&space;\beta_{2,7}\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\beta_{7,1}&space;&&space;\beta_{7,2}&space;&&space;\cdots&space;&&space;\beta_7\\&space;\beta_{8,1}&space;&&space;\beta_{8,2}&space;&&space;\cdots&space;&&space;\beta_{87}\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\beta_{12,1}&space;&&space;\beta_{12,2}&space;&&space;\cdots&space;&&space;\beta_{12,7}\\&space;\end{pmatrix}" title="\begin{pmatrix} \beta_1 & \beta_{1,2} & \cdots & \beta_{1,7}\\ \beta_{2,1} & \beta_2 & \cdots & \beta_{2,7}\\ \vdots & \vdots & \ddots & \vdots \\ \beta_{7,1} & \beta_{7,2} & \cdots & \beta_7\\ \beta_{8,1} & \beta_{8,2} & \cdots & \beta_{87}\\ \vdots & \vdots & \ddots & \vdots \\ \beta_{12,1} & \beta_{12,2} & \cdots & \beta_{12,7}\\ \end{pmatrix}" />

The off-diagonal elements in the first seven rows represents the interactions. The null hypothesis that A vaginae doesn't impact the bias of all the other taxa corresponds to <img src="https://latex.codecogs.com/gif.latex?\beta_{1,2}&space;=&space;\beta_{1,3}&space;=&space;\cdots&space;=&space;\beta_{1,7}" title="\beta_{1,2} = \beta_{1,3} = \cdots = \beta_{1,7}" />, which is a 5 degrees of freedom test. 

We provide two ways to set up this contrasts. The first way is through C.list. Each element in C.list corresponds to one row of the whole contrast matrix (<img src="https://latex.codecogs.com/gif.latex?\mathbb&space;C" title="\mathbb C" /> in the manuscript). Each element in C.list (name it as <img src="https://latex.codecogs.com/gif.latex?C_i" title="C_i" />) is a matrix of the same dimensionality as <img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" />, and puts constraints as <img src="https://latex.codecogs.com/gif.latex?C_i*\beta&space;=&space;0" title="C_i*\beta = 0" />, in which <img src="https://latex.codecogs.com/gif.latex?*" title="*" /> is an element wise product. 

```{r exampleC.list}
C.list = list()
C.list[[1]] = C.list[[2]] = C.list[[3]]  = C.list[[4]] = C.list[[5]] = matrix(0, NROW(mod$beta.matrix), NCOL(mod$beta.matrix))
C.list[[1]][1, 2:3] = c(1, -1)
C.list[[2]][1, c(2, 4)] = c(1, -1)
C.list[[3]][1, c(2, 5)] = c(1, -1)
C.list[[4]][1, c(2, 6)] = c(1, -1)
C.list[[5]][1, c(2, 7)] = c(1, -1)
interTest_Ava =  LLBias.test(mod, C.list = C.list, 
                                n.perm = 1000,  tol = 1e-8)
```

The second approach is through C.matrix, which directly corresponds to <img src="https://latex.codecogs.com/gif.latex?\mathbb&space;C" title="\mathbb C" /> in the manuscript. The number of rows in C.matrix corresponds to the degrees of freedom of the test, and the number of columns of C.matrix is <img src="https://latex.codecogs.com/gif.latex?M&space;\times&space;J" title="M \times J" />. 

```{r}
C.matrix = matrix(0, 5, 84)
C.matrix[1, c(13,25)] = c(1, -1);
C.matrix[2, c(13,37)] = c(1, -1); 
C.matrix[3, c(13,49)] = c(1, -1); 
C.matrix[4,  c(13,61)] = c(1, -1)
C.matrix[5,  c(13,73)] = c(1, -1)
interTest_Ava = LLBias.test(mod, C.matrix = C.matrix, 
                               n.perm = 1000,  tol = 1e-8)
```


## References: 
[1] Zhao N, Satten G. A Log–Linear Model for Inference on Bias in Microbiome Studies. in Press

[2] McLaren MR, Willis AD, Callahan BJ. Consistent and correctable bias in metagenomic sequencing experiments. eLife, 8:e46923, 2019. ISSN 2050-084X.

[3] Brooks JP,Edwards DJ,Harwich MD,Rivera MC,Fettweis JM,Serrano MG,Reris RA,Sheth NU, Huang B, Girerd P, Strauss JF, Jefferson KK, Buck GA. The truth about metagenomics: quantifying and counteracting bias in 16S rRNA studies. BMC Microbiol, 15:66, 2015.

