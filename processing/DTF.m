function [j,DTFX] = DTF(X,f)
%X = EEG_T;
Xf = fft(X);
A = corr(Xf);
H = inv(A);
j = min(abs(Xf-f),2);
%    H{m} = H1(1:m,1:m);1
denom = sum(H(:,1:j(1,1))*H(1:j(1,1),:),2);
DTFX = abs(H(:,j(1,1)))./sqrt(denom);
figure,
pspectrum(DTFX(:)) 
title([{'the frequency '},{int2str(f)},{' Hz is about at time '},{int2str(j(1,1))}, {' ms with DTF approximation'}])
end