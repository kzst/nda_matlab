function [dcorr, dcov]=dcor(varargin)
switch nargin
    case 1 
        X=varargin{1};
        n=size(X,2);
        if isa(varargin{1},'gpuArray')           
            dcorr=gpuArray(zeros(n));
            dcov=gpuArray(zeros(n));
        else
            dcorr=zeros(n);
            dcov=zeros(n);
        end
        for i=drange(1:n)
            for j=drange(i:n)
                [dcorr(i,j),dcov(i,j)]=dCORR(X(:,i),X(:,j));
                dcorr(j,i)=dcorr(i,j);
                dcov(j,i)=dcov(i,j);
            end
        end
    case 2
        P=varargin{1};
        Q=varargin{2};
        if ((min(size(P))==1)&&(min(size(Q))==1))
           [dcorr, dcov] = dCORR(varargin{1},varargin{2}); 
        else
            if (size(P,1)==size(Q,1))
                for i=1:size(Q,2)
                    for j=1:size(P,2)
                        [dcorr(j,i),dcov(j,i)]= dCORR(P(:,j),Q(:,i)); 
                    end
                end
            else
                disp('X and Y must have the same number of rows!');
            end
        end
        
    otherwise
        disp('Two many input paramters!');
end

