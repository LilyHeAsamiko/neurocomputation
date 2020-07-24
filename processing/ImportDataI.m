clear
i = 1;
ampdif = zeros(64, 50); angdif = zeros(64, 50); locdif = zeros(64, 50);
for i = 1:16
    load(strcat('thalamus_source',int2str(i),'.mat'))
    ampdif(i,:) = amplitude_difference;
    angdif(i,:) = angle_difference;
    locdif(i,:) = location_difference;
    clear('amplitude_difference', 'angle_difference', 'location_difference');
end

j = 0;
for i = 17:2:32
        load(strcat('thalamus_source',int2str(i - j),'.mat'))
        ampdif(i,:) = amplitude_difference(1,:);
        ampdif(i+1,:) = amplitude_difference(2,:);
        angdif(i,:) = angle_difference(1,:);
        angdif(i+1,:) = amplitude_difference(2,:);
        locdif(i,:) = location_difference(1,:);
        locdif(i+1,:) = amplitude_difference(2,:);
        j = j + 1;
        clear('amplitude_difference', 'angle_difference', 'location_difference');
end

for i = 33:48
    load(strcat('thalamus_source',int2str(i - 8),'.mat'))
    ampdif(i,:) = amplitude_difference;
    angdif(i,:) = angle_difference;
    locdif(i,:) = location_difference;
    clear('amplitude_difference', 'angle_difference', 'location_difference');
end

j = 0;
for i = 49:2:64
    load(strcat('thalamus_source',int2str(i - 8 - j),'.mat'))
    ampdif(i,:) = amplitude_difference(1,:);
    ampdif(i+1,:) = amplitude_difference(2,:);
    angdif(i,:) = angle_difference(1,:);
    angdif(i+1,:) = amplitude_difference(2,:);
    locdif(i,:) = location_difference(1,:);
    locdif(i+1,:) = amplitude_difference(2,:);
    j = j + 1;
    clear('amplitude_difference', 'angle_difference', 'location_difference');
end

beta = 1.5;
theta0 = [10^-5, 10^-9];
hyperprior = [1, 2];
noise = [0.02, 0.05];
Theta0 = zeros(64,50);
Hyperprior = zeros(64,50);
Beta = repmat(beta, 64,50);
Noise = zeros(64,50);
datatype(1:32,50) = 1;
datatype(33:64,50) = 2;
syn = zeros(64, 50);
for i = 1:32:33
    syn(i:i+7,:) = 1;    
    syn(i+8:i+15,:) = 2;
end    
for i = 1:2:15  
    syn(i+16,:) = 3;
    syn(i+17,:) = 4;
    syn(i+48,:) = 3;
    syn(i+49,:) = 4;
end
for i = 1:8:57
    Theta0(i:i+3,:) = repmat(theta0(1), 4, 50);
    Theta0(i+4:i+7,:) = repmat(theta0(2), 4, 50);
end
for i = 1:4:61
    Hyperprior(i:i+1,:) = repmat(hyperprior(1), 2, 50);
    Hyperprior(i+2:i+3,:) = repmat(hyperprior(2), 2, 50);
end
for i = 1:2:63
    Noise(i,:) = noise(1);
    Noise(i+1,:) = noise(2);
end
% 1. amplitude difference, 2. angle difference, 3. location difference, 4.Beta 
% 5. Hyperprior(1. gamma 2. inverse gamma) 6. Noise(0.02 0.05) 7. localization (1. thalamus 2. somatosensory 3. thalamus_pair 4. somatosensory_pair)
X = [ampdif(:), angdif(:), locdif(:), Beta(:), Hyperprior(:), Noise(:), syn(:)];
Y = Theta0(:);
data = [X,Y];
Y1 = X(:,7);
Y2 = X(:,5);
Y3 = X(:,6);
for i = 1:9
PDCX(i,:) = PDCX(i,:)
end