# nda_matlab
Network-based dimensionality reduction and analysis in MATLAB

This package provides Network-based dimensionality reduction and analysis. 

* Network-based dimensionality reduction and analysis. 

* Dimensional reduction. 
* Plot and biplot functions
* Data generation

#### Author

* Zsolt T. Kosztyan

#### Contributor

* Zsolt T. Kosztyan

#### Maintainer

* Zsolt T. Kosztyan

### Outputs:
L: n by m matrix of factor scores, where n is the number of rows in a datasource, m is tne number of latent factors
C: m by m factor correlation matrix 
COMMUNALITY: n by 1 row vector of communalities
LOADINGS: s by m matrix of factor loadings, where s is the number of selected indicators
LTABLE: s by m table of factor loadings, where s is the number of % selected indicators
S: m by 1 vector of membership

### Input:
data: n by M matrix/table/structure of data source (mandatory)

### Optional input parameters:

XHeader: M by 1 cell array of variable names
CorrMethod|cor_method: Correlation method (optional)
  Pearson|pearson|'1'|1: Pearson's correlation (default)
  Spearman|spearman|'2'|2: Spearman's correlation
  Kendall|kendall|'3'|3: Kendall's correlation
  Distance|distance|'4'|4: Distance correlation
 -otherwise: 1 (Pearson's correlation)
MinCor2|min_R: Minimal square correlation between indicators (default: 0)
MinimalCommunity|min_comm: Minimal number of indicators in a community (default: 2)
Gamma: Gamma parameter in multiresolution null_modell (default: 1)
NullModelType|null_model_type (default: 1);
 NewmannGrivan|'1'|1: Newmann-Grivan's null modell
 AvgDet: Null model is the mean of square correlations between indicators
 MinDet,min_det: Null modell is the specified minimal square correlation
 (min_det)
MinEigCentValue|min_evalue: Minimal EVC value (default: 0.00)
MinCommunality|min_communality: Minimal communality value of indicators (default: 0.25)
ComCommunalities|com_communalities=0.0: Minimal common communalities
RotateMethod: Rotation method (default: none);
Biplots: Draw biplots (default: false)
cuts: Draw correlation graph with cuts value (default: 0 => No
correlation graph)

## Usages:
[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data)
[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data,Xheader)
[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data,Xheader,...)

## Examples:
load CWTS_2020
[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(CWTS_2020)
[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(CWTS_2020,'RotationMethod','varimax','MinimalCommunity',3)

## Requirements:
 Eigenvector centralities (if Matlab release is older than R2020a) 
 (Contributors): 
   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
   Rick Betzel, Indiana University, 2012
   Mika Rubinov, University of Cambridge, 2015

 Modified GenLouvain toolbox (Contributurs):
  Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha, 
 "A generalized Louvain method for community detection implemented in 
 MATLAB," https://github.com/GenLouvain/GenLouvain (2011-2019).
