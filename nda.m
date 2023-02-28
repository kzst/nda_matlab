%Author: Zsolt T. Kosztyan Ph.D habil.
% University of Pannonia, Faculty of Economics, 
%  Department of Quantitative Methods
%----------------
%Implementation of Network-based, non-parametric principal component 
% analysis
%----------------
%Outputs:
%L: n by m matrix of factor scores, where n is the number of rows in a
% datasource, m is tne number of latent factors
%C: m by m factor correlation matrix 
%COMMUNALITY: n by 1 row vector of communalities
%LOADINGS: s by m matrix of factor loadings, where s is the number of 
% selected indicators
%LTABLE: s by m table of factor loadings, where s is the number of 
% selected indicators
%S: m by 1 vector of membership
%---------------- 
%Input:
%data: n by M matrix of data source (mandatory)
%---------------- 
%Optional input parameters:
%---------------- 
%XHeader: M by 1 cell array of variable names
%CorrMethod|cor_method: Correlation method (optional)
%  Pearson|pearson|'1'|1: Pearson's correlation (default)
%  Spearman|spearman|'2'|2: Spearman's correlation
%  Kendall|kendall|'3'|3: Kendall's correlation
%  Distance|distance|'4'|4: Distance correlation
% -otherwise: 1 (Pearson's correlation)
%MinCor2|min_R: Minimal square correlation between indicators (default: 0)
%MinimalCommunity|min_comm: Minimal number of indicators in a community
% (default: 2)
%Gamma: Gamma parameter in multiresolution null_modell (default: 1)
%NullModelType|null_model_type (default: 1);
% NewmannGrivan|'1'|1: Newmann-Grivan's null modell
% AvgDet: Null model is the mean of square correlations between indicators
% MinDet,min_det: Null modell is the specified minimal square correlation
% (min_det)
%MinEigCentValue|min_evalue: Minimal EVC value (default: 0.00)
%MinCommunality|min_communality: Minimal communality value of indicators 
% (default: 0.25)
%ComCommunalities|com_communalities=0.0: Minimal common communalities
%RotateMethod: Rotation method (default: none);
%Biplots: Draw biplots (default: false)
%cuts: Draw correlation graph with cuts value (default: 0 => No
%correlation graph)

%---------------- 
%Usages:
%[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data)
%[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data,Xheader)
%[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data,Xheader,...)
%---------------- 
%Examples:
%load CWTS_2020
%[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(CWTS_2020)
%[L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(CWTS_2020,...
%  'RotationMethod','varimax','MinimalCommunity',3)

%---------------- 

%Requirements:
% Eigenvector centralities (if Matlab release is older than R2020a) 
% (Contributors): 
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012
%   Mika Rubinov, University of Cambridge, 2015

% Modified GenLouvain toolbox (Contributurs):
%  Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha, 
% "A generalized Louvain method for community detection implemented in 
% MATLAB," https://github.com/GenLouvain/GenLouvain (2011-2019).

function [L,C,COMMUNALITY,LOADINGS,LTABLE,S]=nda(data,varargin)
if ispc % In case of Windows OS
    path(path,'.\GenLouvain')  % Modified Louvain's modularity calculations
    addpath(genpath('.\scGEAToolbox-master'))    % Alternative of Louvain's 
    addpath(genpath('.\scGEAToolbox-master\+run'))             % modularity
    if isMATLABReleaseOlderThan('R2022b')
        path(path,'.\EVC')         % Eigen-vector centralities (EVC)
    end
else
    path(path,'./GenLouvain')  % Modified Louvain's modularity calculations
    addpath(genpath('./scGEAToolbox-master'))    % Alternative of Louvain's 
    addpath(genpath('scGEAToolbox-master\+run'))               % modularity
    if isMATLABReleaseOlderThan('R2022b')
        path(path,'./EVC')         % Eigen-vector centralities (EVC)
    end
end % In case of non-Windows OS (Mac, Linux)
tdata=[];
if isstruct(data)
    tdata=struct2table(data);
end
if istable(data)
    tdata=data;
end
if istable(tdata)
    XHeader=tdata.Properties.VariableNames;
    data=double(table2array(tdata));
else
    XHeader=table2cell(table(num2str([1:size(data,2)]')));% Default variable 
                                % names is the count number of the variable
end
data(isnan(data))=0; 
data(isinf(data))=0;
DATA=data;
X=DATA;
cor_method=1; % Deault: Pearson's correlation method
min_R=0;      % Minimal square correlation between indicators (default: 0)
min_comm=2;   % Minimal number of indicators in a community (default: 2)
Gamma=1;      % Gamma parameter in multiresolution null_modell (default: 1)
null_modell_type=1;      % Newmann-Grivan's null modell
min_det=0.25;            % Minimal determination value for specified 
                         %minimal determination null model
min_evalue=0.00;          % Minimal EVC value (default: 0.00)
min_communality=0.00;    % Minimal communality value of indicators 
com_communalities=0.0;   % Minimal common communalities
is_rotated=false;        % No rotation
rotate_method='varimax'; % Default varimax rotation, if rotation is set
biplots=false;           % Default no biplot
cuts=0.0;                % Default no displays of correlation graph
for i=1:nargin-1
    if iscell(varargin{i})
        XHeader=varargin{i}; % Use variable names
        XHeader=strrep(XHeader,'_','\_')
    else
        switch varargin{i}
            case {'CorrMethod','cor_method'} % Correlation method
                switch varargin{i+1}
                    case {'Pearson','pearson','1',1} 
                        cor_method=1;               % Pearson's Correlation
                    case {'Spearman','spearman','2',2}
                        cor_method=2;              % Spearman's Correlation
                    case {'Kendall','kendall','3',3}
                        cor_method=3;               % Kendall's Correlation
                    case {'Distance','distance','4',4}
                        cor_method=4;                % Distance Correlation
                    otherwise
                      disp(['Warning: Only Pearson, Spearman, ',...
                           'Kendall and Distance Correlations are ',...,
                           'implemented. The Pearson correlation will ',...
                           'be set.'])
                        cor_method=1;               % Pearson's Correlation
                end
            case {'MinCor2','min_R'}   % Minimal square correlation between 
                min_R=varargin{i+1};   % indicators 
            case {'MinimalCommunity','min_comm'}        % Minimal number of
                min_comm=varargin{i+1};         % indicators in a community
            case {'Gamma','gamma'}     % Gamma parameter in multiresolution 
                Gamma=varargin{i+1};                          % null_modell
            case {'NullModelType','null_model_type','null_modell_type'}
                switch varargin{i+1}
                    case {'NewmannGrivan','1',1}         % Newmann-Grivan's 
                        null_modell_type=1;                   % null modell
                    case {'AvgDet','2',2}           % Average determination              
                        null_modell_type=2;                   % coefficient
                    case {'MinDet','3',3}           % Minimal determination
                        null_modell_type=3;                   % coefficient
                        min_det=varargin{i+2};
                    otherwise
                        disp('This null model type is not implemented')% If 
                        null_modell_type=1;% the desired null modell is not
                end  % implemented, the default null modell is set to be 1.
            case {'MinEigCentValue','min_evalue'}       % Minimal EVC value
                min_evalue=varargin{i+1};
            case {'MinCommunality','min_communality'} % Minimal communality 
                min_communality=varargin{i+1};        % value of indicators
            case {'ComCommunalities','com_communalitities'}% Minimal common 
                com_communalities=varargin{i+1};            % communalities
            case {'RotateMethod','rotate_method'}           % Rotate method
                is_rotated=true;
                rotate_method=varargin{i+1};
            case {'Biplots','biplots'}      % Draw biplots (default: false)
                biplots=varargin{i+1};
            case {'cuts','disp_r'}  % Draw correlation graph with minimal 
                cuts=varargin{i+1};   % cutting cuts value (default: 0) 
        end
    end
end

switch cor_method
    case 1 % Pearson's correlation
        COR=corr(X);
    case 2 % Spearmen's correlation
        COR=corr(X,'Type','Spearman');
    case 3 % Kendall's correlation
        COR=corr(X,'Type','Kendall');
    case 4 % Distance correlation
        COR=dcor(X);
end

COR(isnan(COR))=0;                       % Correlation matrix of indicators 

R=sparse(min(COR.^2,ones(size(COR))-eye(size(COR))));  % Square correlation
clear COR           % matrix, as an adjecency matrix 

R=(R+R')/2;         % Determination matrix must be symmetric

R=min(R,R>min_R);   % Drop low determinations

% Calculate null modell

kin=sum(R)';
kout=sum(R,2);
l=sum(sum(R));
N=(kout*kin')/l;
N(isnan(N))=0;

% Calculate modularity

coords=ones(size(R,1),1); % Set of indicators (1:included, 0:excluded)

R2=R(coords==1,coords==1);% Part of indicators
MTX=R2;                   % There is no null model
switch null_modell_type  
    case 1                % Null model: multi-resolution version of Newman-
        MTX=R2-N*Gamma;   %-Grivan model  
    case 2                % Null model: Average determination coefficient
        MTX=R2-Gamma*mean(R2(R2>0)); 
    case 3                % Null model: Specified minimal determination
        MTX=R2-Gamma*min_det; %coefficient (default value: 0.5)
        
end

if ispc
    S=louvain(R); % Calculate Louvain's modularity
else
    MTX(MTX<-1)=0; % Dorp low values
    S=genlouvain(MTX,'limit',0); % Calculate Louvain's modularity
end

M=unique(S);

if M(1)==0
    M=M(2:end);
end

for i=1:numel(M)
    if numel(S(S==M(i)))<min_comm
        coords(S==M(i))=0; % Ignore minimal communalities
    end
end

S(coords==0)=0; % Ignore non relevant indicators


n=numel(unique(M)); % Number of communities

Coords=1:numel(S);

if cuts>0
    figure('Name','NDA plot - Modules of correlation graph (without FS)')
    G=graph(R);
    R2=min(R,R>cuts);
    G=graph(R2);
    Coords=1:numel(coords);
    h=plot(G,'k','Layout','force','UseGravity',true);
    labelnode(h,Coords,XHeader);
    HSV=hsv(numel(M));
    for i=1:numel(M)
        highlight(h,Coords(S==M(i)),'NodeColor',HSV(i,:));
        H=R2;
        H(Coords(S~=M(i)),:)=0;
        H(:,Coords(S~=M(i)))=0;
        highlight(h,graph(H),'EdgeColor',HSV(i,:));
    end
end

L=zeros(size(DATA,1),numel(M));
EVCs={};
DATAs={};
for i=1:numel(M)
    Coordsi=Coords((S==M(i))&(coords==1));
    if isMATLABReleaseOlderThan('R2022b')
        EVC=eigenvector_centrality_und(full(R(Coordsi,Coordsi)));
    else
        EVC=centrality(full(R(Coordsi,Coordsi)),"eigenvector");
    end
    EVCs{i}=EVC;
    if numel(EVC(EVC>min_evalue))>2
        L(:,i)=sum(data(:,Coordsi(EVC>min_evalue)).*EVC(EVC>min_evalue)',2);
        coords(Coordsi(EVC<=min_evalue))=0;
        S(Coordsi(EVC<=min_evalue))=0;
    else
        L(:,i)=sum(data(:,Coordsi).*EVC',2);
    end
    EVCs{i}=EVC;
    DATAs{i}=data(:,S==M(i));
end

if ((size(L,2)>1)&&(is_rotated==true))
    L=rotatefactors(zscore(L*pca(L)),'Method','varimax');
else
    L=zscore(L);
end

% Calculate factor loadings and communalities

C=corr(L);
switch cor_method
    case 1 % Pearson's correlation
        LOADINGS=corr(L,data(:,S~=0))';
    case 2 % Spearmen's correlation
        LOADINGS=corr(L,data(:,S~=0),'Type','Spearman')';
    case 3 % Kendall's correlation
        LOADINGS=corr(L,data(:,S~=0),'Type','Kendalln')';
    case 4 % Distance correlation
        LOADINGS=dcor(L,data(:,S~=0))';
end

% Feature selection - using constraint of communality

l=true;    
while l==true
    l=false;
    COMMUNALITY=max(LOADINGS'.^2);
    COMMUNALITY(isnan(COMMUNALITY))=0;
    CoordsS=Coords(S~=0);
    s=S(S~=0);
    coordsS=coords(S~=0);
    
    for i=1:numel(M)
        Coordsi=Coords((S==M(i))&(coords==1));
        CoordsiC=Coords((s==M(i))&coordsS==1);
        COM=COMMUNALITY(CoordsiC);
        com_min=min(COM);
        %%if numel(COM(COM<=min_communality))>2
        if numel(S(S==M(i)))>2
            l=true;
            S(Coordsi(COM<=min_communality))=0;
            coords(Coordsi(COM<=min_communality))=0;
            EVC=EVCs{i};
            EVC=EVC(COM>min_communality);
            EVCs{i}=EVC;
            L(:,i)=sum(data(:,Coordsi(COM>min_communality)).*EVC',2);
        else
            EVC=EVCs{i};
            L(:,i)=sum(data(:,Coordsi).*EVC',2);
        end 
    end
    
    if ((size(L,2)>1)&&(is_rotated==true))
        L=rotatefactors(zscore(L*pca(L)),'Method','varimax');
    else
        L=zscore(L);
    end
    
    % Calculate factor loadings and communalities
    
    C=corr(L);
    switch cor_method
        case 1 % Pearson's correlation
            LOADINGS=corr(L,data(:,S~=0))';
        case 2 % Spearmen's correlation
            LOADINGS=corr(L,data(:,S~=0),'Type','Spearman')';
        case 3 % Kendall's correlation
            LOADINGS=corr(L,data(:,S~=0),'Type','Kendalln')';
        case 4 % Distance correlation
            LOADINGS=dcor(L,data(:,S~=0))';
    end
    COMMUNALITY=max(LOADINGS'.^2);
    if min(COMMUNALITY)>=min_communality
        l=false;
    end
end     


l=false;
while l==false
    l=true;
    CCs=zeros(size(LOADINGS,1),1);
    if size(LOADINGS,2)>1
        CoordsC=Coords(S~=0);
        L2=LOADINGS.^2;
        nL2=size(L2,1);
        for I=1:nL2
            [CJ,~]=max(L2(I,:));
            CJ2=max(L2(I,L2(I,:)~=CJ)); %2nd maximal value;
            if ((CJ>=CJ2+com_communalities)||(CJ>2*CJ2))
            else
                CCs(I)=1; %set of common communalities;
            end
        end
    end
    if numel(CCs(CCs==1))>0 %minimal communality must be eliminated
        Coords_real=CoordsC(CCs==1);
        COM=COMMUNALITY(CCs==1);
        [O_COM,P_COM]=sort(COM);
        Coords_real=Coords_real(P_COM);
        l=true;
        i=1;
        while ((l==true)&&(i<=numel(O_COM)))
            if (numel(S(S==S(Coords_real(i))))>2)
                l=false;
                S(Coords_real(i))=0;
                coords(Coords_real(i))=0;
            end
            i=i+1;
        end
    end
    
    for i=1:numel(M)
        Coordsi=Coords((S==M(i))&(coords==1));
        EVC=eigenvector_centrality_und(full(R(Coordsi,Coordsi)));
        EVCs{i}=EVC;
        L(:,i)=sum(data(:,Coordsi).*EVC',2);
    end
    
    if ((size(L,2)>1)&&(is_rotated==true))
        L=rotatefactors(zscore(L*pca(L)),'Method','varimax');
    else
        L=zscore(L);
    end
    
    % Calculate factor loadings and communalities
    
    C=corr(L);
    switch cor_method
        case 1 % Pearson's correlation
            LOADINGS=corr(L,data(:,S~=0))';
        case 2 % Spearmen's correlation
            LOADINGS=corr(L,data(:,S~=0),'Type','Spearman')';
        case 3 % Kendall's correlation
            LOADINGS=corr(L,data(:,S~=0),'Type','Kendalln')';
        case 4 % Distance correlation
            LOADINGS=dcor(L,data(:,S~=0))';
    end
    COMMUNALITY=max(LOADINGS'.^2);
    
end

if cuts>0
    figure('Name','NDA plot - Modules of correlation graph (with FS)')
    G=graph(R);
    R2=min(R,R>cuts);
    G=graph(R2);
    Coords=1:numel(coords);
    h=plot(G,'k','Layout','force','UseGravity',true);
    labelnode(h,Coords,XHeader);
    HSV=hsv(numel(M));
    for i=1:numel(M)
        highlight(h,Coords(S==M(i)),'NodeColor',HSV(i,:));
        H=R2;
        H(Coords(S~=M(i)),:)=0;
        H(:,Coords(S~=M(i)))=0;
        highlight(h,graph(H),'EdgeColor',HSV(i,:));
    end
end


NDA=LOADINGS;                        % Structured output of factor loadings
LTABLE=array2table(NDA,'RowNames',XHeader(S~=0));

if biplots==true
    if n>1
        figure('Name','2D biplots')
        k=1;
        for i=1:n
            for j=1:n
                subplot(n,n,k)
                if (i==j)
                    histfit(LOADINGS(:,i),10,'kernel')
                    title(strcat('NDA',num2str(i)));                    
                else
                    biplot(LOADINGS(:,[i,j]),'Scores',L(:,[i,j]),'VarLabels',XHeader(S~=0))
                    xlabel(['Latent variable ',num2str(i)]); ylabel(['Latent variable ',num2str(j)]);
                end
                k=k+1;
            end
        end
        if n>2
            figure('Name','3D biplots')
            for i=3:n
                subplot(1,n-2,i-2);
                biplot(LOADINGS(:,i-2:i),'Scores',L(:,i-2:i),'VarLabels',XHeader(S~=0))
                xlabel(['Latent variable ',num2str(i-2)]); ylabel(['Latent variable ',num2str(i-1)]);zlabel(['Latent variable ',num2str(i)]);
            end
        end
    end
end

