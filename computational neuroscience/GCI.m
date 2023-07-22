function [GCIM,pGCIM] = GCI(xM,m,maketest)
% [GCIM,pGCIM] = GCI(xM,m,maketest)
% GCI computes the Granger Causality index (GCI) for all time series 
% pairs (Xi,Xj) in the vector time series given in 'xM', for both 
% directions Xi->Xj and Xj->Xi. For the computation of GCI the unrestricted 
% and restricted VAR and AR models of the given order 'm' are fitted for 
% each response variable. Further, if it is selected (maketest=1) the 
% parametric significance test is performed for each pair of variables.
% INPUTS
% - xM          : the vector time series of size n x K 
% - m           : the order of the restricted AR and unrestricted VAR model 
% - maketest    : If 1 make parametric significance test and give out 
%                 the p-value for each pair of variables
% OUTPUTS
% - GCIM       : The matrix of size K x K of the GCI values. Each cell
%                (i,j) has the value GCI(Xi->Xj). 
% - pGCIM      : The matrix of size K x K of the p-values of the 
%                parametric significance test for GCI (using the F-statistic, see Econometric 
%                Analysis, Greene, 7th Edition, Sec 5.5.2)
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

% Fit the AR(m) model for each response and compute the sum of squared 
% errors (SSE) in order to use it later as the SSE for the restricted 
% model.
xrss0V = NaN(K,1);
for iK=1:K
    xiM = xxM(:,[(iK-1)*m+1:iK*m]);
    xbV = xiM\xM(m+1:n,iK);
    xpreV = xiM*xbV;
    xerrV = xM(m+1:n,iK) - xpreV;
    xrss0V(iK) = sum(xerrV.^2);
end

GCIM = NaN(K,K);
pGCIM = NaN(K,K);
for iK=1:K
    for jK=iK+1:K
        xiM = xxM(:,[(iK-1)*m+1:iK*m]);
        xjM = xxM(:,[(jK-1)*m+1:jK*m]);        
        % Full regression on Xj from Xi-past and Xj-past (unrestricted)
        xbV = [xiM xjM]\xM(m+1:n,jK);
        xpreV = [xiM xjM]*xbV;
        xerrV = xM(m+1:n,jK) - xpreV;
        xjrss = sum(xerrV.^2);
        GCIM(iK,jK) = log(xrss0V(jK)/xjrss);
        % Full regression on Xi from Xi-past and Xj-past (unrestricted)
        xbV = [xiM xjM]\xM(m+1:n,iK);
        xpreV = [xiM xjM]*xbV;
        xerrV = xM(m+1:n,iK) - xpreV;
        xirss = sum(xerrV.^2);
        GCIM(jK,iK) = log(xrss0V(iK)/xirss);
        if maketest
            ndf = n-m-2*m;
            xfstat = ((xrss0V(jK)-xjrss)/m)/(xjrss/ndf); 
            pGCIM(iK,jK) = 1 - fcdf(xfstat,m,ndf);
            xfstat = ((xrss0V(iK)-xirss)/m)/(xirss/ndf); 
            pGCIM(jK,iK) = 1 - fcdf(xfstat,m,ndf);
        end
    end
end
