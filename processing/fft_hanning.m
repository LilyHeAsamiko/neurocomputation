function FFT_M = fft_hanning(X,M,fs,fpass,fstop,fc)
%X = zef_EEG.measurements(33,:);
X = fft(EEG_T);
%pspectrum(X);
n = 2*M+1;
w = [2*fc*sinc(2*pi*(1:n)) zeros(1,N-n)]; %hanning window
X_W = zeros(72,N+1); %initialte the zero pad out
FFT_M = zeros(72,360);
X_W(1,:) = [X(1,1:length(w)+1)] .* fft([2*fc,w]);% fft with convolution   
M_X(1,:) = mean(X_W(1,:),2); %average the slide windows
FFT_M(1,1) = mean(M_X(1,:));
for i = 2:size(X,1)
    for j = 1:360-N
    X_W(i,j:j+N) = X(i,1:length(w)+1) .* fft([2*fc,w]);% fft with convolution   
    M_X(i,j) = mean(X_W(i,j:j+N),2); %average the slide windows
    
    end
    FFT_M(i,j) = mean(M_X(i,:));
    FFT_M(:,j) = (FFT_M(:,j)+FFT_M(:,j+1))/2;
end

% Display real part of windowed signal and the Hanning window:
subplot(3,1,1)
plot(1:N,w(1,:),'-');
title('Hanning Window and Windowed, Zero-Padded');
xlabel('Time (samples)'); ylabel('Amplitude');

subplot(3,1,2)
plot(1:360,real(FFT_M(1,:)));
title('Fast Fourier (Real Part)');
xlabel('Time (samples)'); ylabel('frequency(Hz)');

subplot(3,1,3)
plot(1:360,(abs(FFT_M(1,:))).^2,'r-');
title('Fast Fourier (Real Part)');
xlabel('Time (samples)'); ylabel('Spectrum');
hold on
plot(1:360,(fft(X(1,:))).^2,'b-')
end