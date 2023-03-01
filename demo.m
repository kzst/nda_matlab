%% Demonstration of Netrwork-based Dimensionality Reduction and Analysis
% (NDA) method

%Author: Zsolt T. Kosztyan Ph.D habil.
% University of Pannonia, Faculty of Economics, 
%  Department of Quantitative Methods

%% Simple dimension reduction without Feature Selection

load CWTS_2020 % Loading CWTS LEIDEN RANKING DATA (1176 HEI x 42 variable)

[Scores,CMTX,COMMUNALITY,LOADINGS,LTABLE,MEMBERSHIPS]=nda(CWTS_2020);

disp('Print Loading Table')

LTABLE

%% Using Spearman's Correlation, Plot Biplot

%cor_method=2 %Employing Spearman's corelation
%biplots=true %Plot Biplots

[Scores,CMTX,COMMUNALITY,LOADINGS,LTABLE,MEMBERSHIPS]=nda(CWTS_2020,...
    'cor_method',2,'biplots',true);

disp('Print Loading Table')

LTABLE

%% Using Spearman's Correlation with FS, Plot Network Graph

%cor_method=2 %Employing Spearman's corelation
%cuts=0.2 %Plotting correlation graph, with minimal edge 0.2 value
%min_evalue=0.2 %Setting minimal eignevector centrality value to 0.2

[Scores,CMTX,COMMUNALITY,LOADINGS,LTABLE,MEMBERSHIPS]=nda(CWTS_2020,...
    'cor_method',2,'cuts',0.2,'min_evalue',0.2);

disp('Print Loading Table')

LTABLE