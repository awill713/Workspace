
clear
seed = 5;
rng(seed)

sampleRate = 400e3;

%Design stimulus parameters
freqLow = 440;
freqHigh = 880;
uniqueTones = 10;

samRate = 1; %Hz
samMagnitude = [0 0.5 1]; %max percent reduction in amplitude

repeats = 10;

toneDuration = 5000; %milliseconds
ITI = 1000; %milliseconds
rampDuration = 5; %milliseconds

%Load filter
% filterName = 'E:\calibration\Filters\20180427_2Pspkr_IC_NidaqInvFilt_3k-80k_fs400k.mat';
% load(filterName);

%Create stimulus sequence
logToneRange = linspace(log10(freqLow),log10(freqHigh),uniqueTones);

index = zeros(uniqueTones*length(samMagnitude),3);
for i = 1:size(index,1)
    index(i,1) = i;
    index(i,2) = 10^logToneRange(floor((i-1)/3)+1);
    index(i,3) = samMagnitude(mod(i-1,3)+1);
end

toneOrder = [];
for r = 1:repeats
    trialOrder = datasample(index(:,1),size(index,1),'Replace',false)';
    toneOrder = [toneOrder trialOrder];
end

%Create stimulus vector
toneSamples = toneDuration*sampleRate/1000;
ITIsamples = ITI*sampleRate/1000;
presentationSamples = toneSamples + ITIsamples;
rampSamples = rampDuration * sampleRate/1000;

stimulus = zeros(1,length(toneOrder)*presentationSamples);

for i = 1:length(toneOrder)
        tempFreq = index(toneOrder(i),2);
        coef = tempFreq*2*pi/sampleRate;
        
        tempSignal(1:toneSamples) = cos(coef*(1:toneSamples));
        
        tempSignal = applyRamp_AMW(tempSignal,rampSamples);
        
        samCoef = samRate*2*pi/sampleRate;
        m = 0.5*index(toneOrder(i),3);
        b = 1-m;
        samSignal(1:toneSamples) = m*cos(samCoef*(1:toneSamples)) + b;

        tempSignal = tempSignal .* samSignal;
        
        sampleIndex = (i-1)*presentationSamples+1;
        stimulus(1,sampleIndex:(sampleIndex+toneSamples-1)) = tempSignal;
%         stimulus(2,sampleIndex:(sampleIndex+toneSamples-1)) = 5;
        i
end

% stimulus(1,:) = conv(stimulus,'FILT','same');

% stimulus = stimulus';

minutes = floor(length(stimulus)/sampleRate/60);
sec = (length(stimulus)/sampleRate/60 - minutes)*60;
disp(['Stimulus is ' num2str(minutes) ' minutes and ' num2str(sec) ' seconds']);

%Save files
audiowrite(['SAMtest.wav'],stimulus,sampleRate);
% 
% stimInfo.fileName = fileName;
% stimInfo.sampleRate = sampleRate;
% stimInfo.frequencies = 10.^logToneRange;
% stimInfo.uniqueTones = uniqueTones;
% stimInfo.samMagnitudes = samMagnitudes;
% stimInfo.samRate = samRate;
% stimInfo.toneDuration = toneDuration; % tone duration in ms
% stimInfo.ITI = ITI; % inter tone interval duration in ms
% stimInfo.repeats = repeats; % number of repeats of each tone
% stimInfo.filterName = filterName;
% stimInfo.order = toneOrder;
% stimInfo.stimGenFunc = 'toneGenerator.m';
% 
% save([fileName '_stimInfo.mat'],'stimInfo')

