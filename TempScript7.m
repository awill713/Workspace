
sampleRate = 40e3;

toneDuration = 5000;
toneSamples = toneDuration*sampleRate/1000;

fundFreq = 440;
freqMax = 20000;
numFreq = 10;
overtoneAttn = 0.5;

logRamp = logspace(-1,0,toneSamples);
linRamp = linspace(0.5,1.5,toneSamples);

freqList = logspace(log10(fundFreq),log10(freqMax),numFreq);
freqList = fundFreq:fundFreq:freqMax;

stimulus = zeros(1,toneSamples);

for f=1:length(freqList)
    coef = freqList(f)*2*pi/sampleRate;
    period = sampleRate/freqList(f);
    phaseOffset = randi(round(period),1);
    tempStim = cos(coef*((1+phaseOffset):(toneSamples+phaseOffset))) * overtoneAttn^(f-1);
%     tempStim = cos(coef*(1:toneSamples)) * overtoneAttn^(f-1);
    stimulus = stimulus + tempStim;
end
stimulus = stimulus./max(abs(stimulus));
stimulus = stimulus.*logRamp;

figure;
spectrogram(stimulus,[],[],[],sampleRate,'yaxis');
ylim([0 20]);

figure;periodogram(stimulus,[],[],sampleRate);

figure;plot(stimulus);

sound(stimulus,sampleRate);

%%

%Design stimulus parameters
freq1 = 440;
freq2 = 880;
freq3 = 660;

toneDuration = 5000; %milliseconds

toneSamples = toneDuration*sampleRate/1000;

logramp = logspace(-1,0,toneSamples);
linramp = linspace(0,1,toneSamples);

stim1 = zeros(1,toneSamples);
stim2 = zeros(1,toneSamples);
stim3 = zeros(1,toneSamples);
stim4 = zeros(1,toneSamples);

coef1 = freq1*2*pi/sampleRate;
coef2 = freq2*2*pi/sampleRate;
coef3 = freq3*2*pi/sampleRate;

stim1(1:toneSamples) = cos(coef1*(1:toneSamples));
stim2(1:toneSamples) = cos(coef2*(1:toneSamples));
stim3(1:toneSamples) = cos(coef3*(1:toneSamples));
stim4(1:toneSamples) = rand(1,toneSamples);

stim = (stim1+stim2+stim3)/3;
sound(stim.*logramp,sampleRate)
% sound(stim.*linramp,sampleRate)

%%
saw = sawtooth(0:0.1:1000*pi);
figure;plot(saw)
figure;
periodogram(saw)