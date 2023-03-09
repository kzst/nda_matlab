function [ dCorr, dCov] = dCORR( varargin )
%Calculation of the Improved version of distance correlation
%Created by: Zsolt T. Kosztyán, Attila I. Katona, University of Pannonia,Department of Quantitative Methods, 2014.
%--------------
% REFERENCES:
%--------------
% Székely G., Rizzo M. (2009) "Brownian distance covariance", The annals of Applied Statistics
% Székely G., Rizzo M. (2013) "The distance correlation t-test of independence in high dimension", JMVA
%--------------
%  INPUTS:
%--------------
% x and y are the input vectors
%--------------
%OUTPUTS:
%--------------
%1. dCorr: Distance correlation of x and y , (-1 <= dCorr <= 1)
%2. dCov: Distance covariance of x and y
%--------------
%USAGE
%--------------
% set the path of the dCORR function as current folder
% type:  [ dCov, dVarx, dVary, dCorr ] = dCORR( x,y )

switch nargin
    case 1
        XY=varargin{1};
        if isa(XY,'table');
            XY=table2array(XY);
        end
        if isa(XY,'cell');
            XY=cell2mat(XY);
        end
        xy=reshape(XY,[],2);
        x=xy(:,1);
        y=xy(:,2);

    case 2
        x=varargin{1};
        y=varargin{2};
end

n=numel(x);% calculating the length of x
if isa(x,'gpuArray')
    A=zeros(n,'gpuArray');
    A2=zeros(n,'gpuArray');
    B=zeros(n,'gpuArray');
    B2=zeros(n,'gpuArray');
    a=zeros(n,'gpuArray');
    b=zeros(n,'gpuArray');
else
    A=zeros(n);
    A2=zeros(n);
    B=zeros(n);
    B2=zeros(n);
    a=zeros(n);
    b=zeros(n);
end
%-------------------------------------------------------
%Check if the sizes of the inputs match
%-------------------------------------------------------
if numel(x) ~= numel(y)
    error('Inputs must have the same number of rows') %gives an error message
end
%----------------------------------------------------------------
% Calculate doubly centered distance matrix for x (Article 2009)
%----------------------------------------------------------------
for i=1:n
    for j=1:n
        a(i,j)=pdist2(x(i),x(j),'euclidean'); % calculating euclidean distance between each point in x
    end
end
colmean_a=mean(a,1); % column mean of the distance matrix a
rowmean_a=mean(a,2); % row mean of the distance matrix a
mean_a=mean2(a)*ones(numel(x),1); % the mean vector of the distance matrix a

for i=1:n
    for j=1:n
        A(i,j)=a(i,j)-colmean_a(i)-rowmean_a(j)+mean_a(i);
    end
end
%----------------------------------------------------------------
% Calculate doubly centered distance matrix for y (Article 2009)
%----------------------------------------------------------------
for i=1:n
    for j=1:n
        b(i,j)=pdist2(y(i),y(j),'euclidean'); % calculating euclidean distance between each point in y
    end
end
colmean_b=mean(b,1); % column mean of the distance matrix b
rowmean_b=mean(b,2); % row mean of the distance matrix b
mean_b=mean2(b)*ones(numel(y),1); % the mean vector of the distance matrix b

for i=1:n
    for j=1:n
        B(i,j)=b(i,j)-colmean_b(i)-rowmean_b(j)+mean_b(i);
    end
end
%-------------------------------------------------------
% Correction of A and B (Article 2013)
%-------------------------------------------------------
for i=1:n
    for j=1:n
        if (i==j)
            A2(i,j)=n/(n-1)*(colmean_a(i)-mean_a(i)); % correction of the diagonal vector in A
        else
            A2(i,j)=n/(n-1)*(A(i,j)-(a(i,j)/numel(x)));
        end
    end
end

for i=1:n
    for j=1:n
        if (i==j)
            B2(i,j)=n/(n-1)*(colmean_b(i)-mean_b(i)); % correction of the diagonal vector in B
        else
            B2(i,j)=n/(n-1)*(B(i,j)-(b(i,j)/n));
        end
    end
end
%-----------------------------------------------------------------
%Calculation of distance covariance, and variances (Article 2013)
%-----------------------------------------------------------------

%Distance covariance
ABcov=A2.*B2;
diagABcov=diag(ABcov);
ABcov(logical(eye(size(ABcov)))) = 0; %set the values of the diagonal vector to zero
Ucov=sum(sum(ABcov))-(2/(n-2))*sum(diagABcov);
dCov=Ucov/((n*(n-3)));
%Distance variance for x
AAvarx=A2.*A2;
diagAAvarx=diag(AAvarx);
AAvarx(logical(eye(size(AAvarx)))) = 0;%set the values of the diagonal vector to zero
Uvarx=sum(sum(AAvarx))-(2/(n-2))*sum(diagAAvarx);
dVarx=Uvarx/((n*(n-3)));
%Distance variance for y
BBvary=B2.*B2;
diagBBvary=diag(BBvary);
BBvary(logical(eye(size(BBvary)))) = 0;%set the values of the diagonal vector to zero
Uvary=sum(sum(BBvary))-(2/(n-2))*sum(diagBBvary);
dVary=Uvary/((n*(n-3)));
%----------------------------------------------------
%Calculation of Distance Correlation
%----------------------------------------------------
if dVarx*dVary>0
    dCorr=dCov/(sqrt(dVarx*dVary));
else
    dCorr=0;
end
end