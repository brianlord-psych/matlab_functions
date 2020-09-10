% Example implementation of complexity function
% EEGlab is required to be open

load sampleEEGdata;
data = EEG.data;

n = 64; % number of electrodes
T = 10; % number of trials to average

% Initialize matrices
Fractals = zeros(n,T);
SampleEntropy = zeros(n,T);
SpectralEntropy = zeros(n,T);

HF = zeros(1,n);
SaE = zeros(1,n);
SpE = zeros(1,n);

% Calculate complexity measures for all conditions
for a = 1:n
    for t = 1:T
        Fractals(a,t) = complexity(data(a,:,t),'HFD',5);
        SampleEntropy(a,t) = complexity(data(a,:,t),'SE',3,0.2);
        SpectralEntropy(a,t) = complexity(data(a,:,t),'HHSE');
    end
    
    % Take mean across trials
    HF(a) = mean(Fractals(a,:))
    SaE(a) = mean(SampleEntropy(a,:));
    SpE(a) = mean(SpectralEntropy(a,:));
end

% Normalization of measures for comparison between them
norm_comp = zeros(64,3);

norm_comp(:,1) = (HF-mean(HF))/std(HF);
norm_comp(:,2) = (SaE-mean(SaE))/std(SaE);
norm_comp(:,3) = (SpE-mean(SpE))/std(SpE);

% Plot
figure(11);

subplot(411);
plot(SaE(:));
title('Sample Entropy');
xlabel('EEG electrodes');
xlim([1 64]);

subplot(412);
plot(SpE(:));
title('Hilbert-Huang Spectral Entropy');
xlabel('EEG electrodes');
xlim([1 64]);

subplot(413);
plot(HF(:));
title('Higuchi Fractal Dimension');
xlabel('EEG electrodes');
xlim([1 64]);

subplot(414);
plot(norm_comp);
title('Normalized Comparison of All Measures');
xlabel('EEG electrodes');
xlim([1 64]);

% Topoplot
figure(12);
subplot(221);
topoplot(norm_comp(:,3),EEG.chanlocs,'maplimits','minmax');
title('Sample Entropy');

subplot(222);
topoplot(norm_comp(:,2),EEG.chanlocs,'maplimits','minmax');
title('Hilbert-Huang Spectral Entropy');

subplot(223);;
topoplot(norm_comp(:,1),EEG.chanlocs,'maplimits','minmax');
title('Higuchi Fractal Dimension');




% Running Window
N = length(data);
W = 300;

for x = 1:N-W
    compF(x) = complexity(data(1,x:(x+W)),'HFD');
    compS(x) = complexity(data(1,x:(x+W)),'SE');
    compH(x) = complexity(data(1,x:(x+W)),'HHSE');
end

compF = (compF-mean(compF))/std(compF);
compS = (compS-mean(compS))/std(compS);
compH = (compH-mean(compH))/std(compH);
figure(16); hold on;
plot(compF);
plot(compS);
plot(compH);

%%%
figure(13);
for a=1:7
    subplot(7,1,a);
    plot(imf(:,a));
end