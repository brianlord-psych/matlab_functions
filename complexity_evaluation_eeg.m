signal = EEG.data(1,:,1);

%% Complexity calcs for EEG signal
% HFD specs
kmax = 384;
HFDvalues = NaN(1,kmax-1);

% SampEn specs
dimmax = 50;
r = 0.15;
SampEnvalues = NaN(1,dimmax-1);


% struct for all values
values = struct('fx',signal, 'SampEn',NaN(length(dimmax)), 'HFD',NaN(length(kmax)));

% kmax effects
for k = 2:kmax
    HFDvalues(k) = complexity(signal,'HFD',k);
    disp(k);
end

% dim effects
for d = 2:dimmax
    SampEnvalues(d) = FSampEn(signal,r,d);
    disp(d);
end

values.SampEn = SampEnvalues;
values.HFD = HFDvalues;

values.SampEn(values.SampEn>2) = NaN;

%% Plots

% SampEn plot
figure(1)
t1=tiledlayout(2,1)
nexttile
plot(values.SampEn)
ylabel('SampEn')
xlabel('Embedding dimension')
nexttile
plot(signal);



% HFD Plot
figure(2)
t2=tiledlayout(2,1)
nexttile
plot(values.HFD)
ylabel('HFD')
xlabel('kmax')
nexttile
plot(signal)
title(t2,'Higuchi Fractal Dimension');


