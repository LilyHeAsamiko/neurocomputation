function FFT_M = fft_slide(X,M1,fs,fpass,fstop,fc,t1,t2)
%X = zef_EEG.measurements(33,:);
w1 = zeros(1,t2-t1);
w1 = t1:t2; %slide window 
N1 = t2-t1;
X_W1(t1:t2) = fft(X(t1:t2)).* fft(w1);
M_X1(t1) = mean(X_W1(t1:t2),2);
FFT_M(1) = mean(M_X1(t1:t1+1));
for j = 2:360-N1
    X_W1(j:j+N1) = fft(X(j:j+length(w1)-1)) .* fft(w1);% fft with convolution   
    M_X1(j) = mean(X_W1(j:j+N1),2); %average the slide windows
    FFT_M(j) = fft(M_X1(j));
end
% Display real part of windowed signal and the Hanning window:

figure,
subplot(2,1,1)
plot(1:N1+1,w1(:),'*');
title('Slide Window(from t1 to t2)');
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2)
plot(1:360,real(fft(X(:,31:390))).*imag(fft(X(:,31:390))),'b-');
plot on
plot((1:360-N1),log*(real(FFT_M(:)).*imag(FFT_M(:))),'r-');
title('Fast Fourier (Power Spectrum )with slide window (from t1 to t2)');
hold off

figure,
plot(1:360,(abs(fft(X(:,31:390)))).^2,'b-');
plot on
plot((1:360-N1),log((abs(FFT_M(:))).^2),'r-');
title('Fast Fourier log(Power Spectrum )with slide window (from t1 to t2)');

end