% cross spectrum matrix computation
function [CS,E,SNR] = cross_spectral(X,f,e)
%X = zef.measurements(e,:);
Xf = fft(X);
j = min(abs(Xf-f),2);
Af= corr(Xf(:, j(1,1)));
CS = Xf*Xf'./cov(X(e,:));
showm = CS(:);
ind = 1:size(CS,2):size(CS,2)^2;
E = Xf.*Af;
SNR = 10*log(CS/(E*E'));
figure,
spectrogram(showm(ind));
title('(Sensor) Cross sprctral');
end

