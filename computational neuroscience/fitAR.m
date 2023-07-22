function resV = fitAR(xV,p)
% resV = fitAR(xV,p)
% It gives the array of residuals from the fitting of AR(p) on the time
% series xV. The AR(p) is fitted without using the ident toolbox
epsilon = 10^(-7);
n = length(xV);
mx = mean(xV);
xxV = xV-mx;
xM = NaN(n-p,p);
for ip=1:p
    xM(:,ip) = xxV(p-ip+1:n-ip);
end
yV = xxV(p+1:n);
[Ux, Sx, Vx] = svd(xM, 0);
q = length(find(diag(Sx)>epsilon));
tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
lsbV = tmpM * yV;
xpreV = xM * lsbV;
xpreV = [xxV(1:p);xpreV];
resV = xxV - xpreV;
