function [CGCIM,pCGCIM] = CGCI(xM,m,maketest)
% [CGCIM,pCGCIM] = CGCI(xM,m,maketest)
% CGCI computes the conditional Granger Causality index (CGCI) for all 
% time series pairs (Xi,Xj) in the presence of the rest time series in the
% vector time series given in 'xM', for both directions (Xi->Xj|(X-{Xi,Xj} 
% and Xj->Xi|(X-{Xi,Xj}). For the computation of CGCI the unrestricted and
% restricted VAR models of the given order 'm' are fitted for each 
% response variable. Further, if it is selected (maketest=1) the parametric
% significance test is performed for each pair of variables.
% INPUTS
% - xM          : the vector time series of size n x K 
% - m           : the order of the restricted and unrestricted VAR model 
% - maketest    : If 1 make parametric significance test and give out 
%                 the p-value for each pair of variables
% OUTPUTS
% - CGCIM       : The matrix of size K x K of the CGCI values. Each cell
%                 (i,j) has the value CGCI(Xi->Xj). 
% - pCGCIM      : The matrix of size K x K of the p-values of the 
%                 parametric significance test for CGCI (using the F-statistic, see Econometric 
%                 Analysis, Greene, 7th Edition, Sec 5.5.2)
if nargin==2
    maketest = 0;
end
[n,K] = size(xM);
xM = xM - repmat(mean(xM),n,1);

% The lag matrix (at each row) for system X
xxM = NaN*ones(n-m,K*m);
for iK=1:K
    colnow = (iK-1)*m+1;
    for im=1:m
        xxM(:,colnow+m-im) = xM(im:n-1-m+im,iK);
    end
end

xrss1V = NaN*ones(K,1);
for iK=1:K
    xbV = xxM\xM(m+1:n,iK);
    xpreV = xxM*xbV;
    xerrV = xM(m+1:n,iK) - xpreV;
    xrss1V(iK) = sum(xerrV.^2);
end
CGCIM = NaN*ones(K,K);
pCGCIM = NaN*ones(K,K);
for iK=1:K
    for jK=iK+1:K
        % Autoregression on Xj from all but Xi
        irestV = setdiff([1:K*m],[(iK-1)*m+1:iK*m]);
        xbV = xxM(:,irestV)\xM(m+1:n,jK);
        xpreV = xxM(:,irestV)*xbV;
        xerrV = xM(m+1:n,jK) - xpreV;
        xrss0 = sum(xerrV.^2);
        CGCIM(iK,jK) = log(xrss0/xrss1V(jK));
        % Autoregression on Xi from all but Xj
        jrestV = setdiff([1:K*m],[(jK-1)*m+1:jK*m]);
        ybV = xxM(:,jrestV)\xM(m+1:n,iK);
        ypreV = xxM(:,jrestV)*ybV;
        yerrV = xM(m+1:n,iK) - ypreV;
        yrss0 = sum(yerrV.^2);
        CGCIM(jK,iK) = log(yrss0/xrss1V(iK));
        % Compute p-values of the parametric F-test
        if maketest
            ndf = n-m-K*m; 
            xfstat = ((xrss0-xrss1V(jK))/m)/(xrss1V(jK)/ndf); 
            pCGCIM(iK,jK) = 1 - fcdf(xfstat,m,ndf);
            yfstat = ((yrss0-xrss1V(iK))/m)/(xrss1V(iK)/ndf); 
            pCGCIM(jK,iK) = 1 - fcdf(yfstat,m,ndf);
        end
    end
end
