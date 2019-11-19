
% This script makes a stimulus of pure tones of various amplitudes with or
% without opto laser (simultaneous onset). Shorter interburst
% interval because it's designed for ephys, not 2-photon.

clear
seed = 5;
rng(seed)

fileName = ['D:\stimuli\Aaron\Opto\' datestr(now,'yyyymmdd') '_multiAmpTones_opto_1msPulse_ONE_REPEAT_70dB_400k_' sprintf('%03d',seed) ]; %don't append .wav or .mat

sampleRate = 400e3;

toneDuration = 100; %ms
optoPulseDuration = 1; %ms
repeats = 1;
IBI = 500; %interburst interval, ms
freqLow= 4e3;
freqHigh = 70e3;
uniqueTones = 20;
intensities = [40 50 60 70 80]; %db
rampDuration = 5; %ms
preStimSilence = 1000; %ms

%Load filter
filterName = 'D:\GitHub\filters\20191029_ephysBoothSpkrRight_300-80k_fs400k.mat';
load(filterName);

%Prepare stimulus sequence
logToneRange = linspace(log10(freqLow),log10(freqHigh),uniqueTones);
uniqueEvents = uniqueTones*length(intensities)*2; %laser on and laser off
index = zeros(uniqueEvents,4); %index,frequency,intensity,laser
for i = 1:size(index,1)
    index(i,1) = i; %index
    index(i,2) = 10^logToneRange(mod(i-1,uniqueTones)+1); %frequency
    index(i,3) = intensities(floor((mod(i-1,uniqueEvents/2))/uniqueTones+1)); %intensity
    index(i,4) = floor(i/(uniqueEvents/2)); %laser
end
order = repmat(index(:,1)',[1 repeats]);
order=order(randperm(length(order)));

%Create stimulus
totalSamples = (repeats*uniqueEvents*(toneDuration+IBI)+preStimSilence)*sampleRate/1000;
stimulus = zeros(3,totalSamples);

toneSamples = toneDuration*sampleRate/1000;
IBIsamples = IBI*sampleRate/1000;
preStimSilenceSamples = preStimSilence*sampleRate/1000;
rampSamples = rampDuration*sampleRate/1000;
pulseSamples = optoPulseDuration*sampleRate/1000;

for i = 1:length(order)
    freq = index(order(i),2);
    volume = index(order(i),3);
    laserStatus = index(order(i),4);
    
    coef = freq*2*pi/sampleRate;
    tempSignal(1:toneSamples) = cos(coef*(1:toneSamples));
    tempSignal = applyRamp_AMW(tempSignal,rampSamples);
    
    attn = 70-volume;
    tempSignal = tempSignal.*10^(-attn/20);
    
    pulse = zeros(1,pulseSamples*2);
    pulse(1:pulseSamples) = 5;
    pulseCount = ceil(toneSamples/pulseSamples);
    pulse = repmat(pulse,[1 pulseCount]);
    pulse=pulse(1:toneSamples);

    sampleIndex = (i-1)*(toneSamples+IBIsamples)+1+preStimSilenceSamples;
    stimulus(1,sampleIndex:(sampleIndex+toneSamples-1)) = tempSignal;
    stimulus(2,sampleIndex:(sampleIndex+toneSamples-1)) = 5;

    if laserStatus %laser present
        stimulus(3,sampleIndex:(sampleIndex+toneSamples-1)) = pulse;
    end
end

stimulus(1,:) = conv(stimulus(1,:),FILT,'same');

stimulus = stimulus'./10;

minutes = floor(length(stimulus)/sampleRate/60);
sec = (length(stimulus)/sampleRate/60 - minutes)*60;
disp(['Stimulus is ' num2str(minutes) ' minutes and ' num2str(sec) ' seconds']);

%Save files
audiowrite([fileName '.wav'],(stimulus),sampleRate);


stimInfo.seed = seed;
stimInfo.fileName = fileName;
stimInfo.sampleRate = sampleRate;
stimInfo.toneDuration = toneDuration; % noise duration in ms
stimInfo.optoPulseDuration = optoPulseDuration;
stimInfo.frequencies = 10.^logToneRange;
stimInfo.uniqueTones = uniqueTones;
stimInfo.repeats = repeats;
stimInfo.index = index;
stimInfo.order = order;
stimInfo.preStimSilence = preStimSilence;
stimInfo.intensities = intensities; % in dB
stimInfo.IBI = IBI; % inter stim interval duration in ms
stimInfo.filter = filterName;

save([fileName '_stimInfo.mat'],'stimInfo');