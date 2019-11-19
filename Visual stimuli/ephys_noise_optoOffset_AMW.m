
% This script makes a stimulus of noise bursts (not clicks) with opto laser
% at various timing offsets relative to the sound onset. Shorter interburst
% interval because it's designed for ephys, not 2-photon.

clear
seed = 5;
rng(seed)

fileName = ['D:\stimuli\Aaron\Opto\' datestr(now,'yyyymmdd') '_noise_optoOffset_25msPulse_50rep_70dB_400k_' sprintf('%03d',seed) ]; %don't append .wav or .mat

sampleRate = 400e3;

burstDuration = 100; %ms
optoOffset = [-20 -10 -5 0 5 10 20]; %ms
optoPulseDuration = 25; %ms
repeats = 50;
IBI = 500; %interburst interval, ms
intensity = 70; %db
rampDuration = 5; %ms
preStimSilence = 1000; %ms

%Load filter
filterName = 'D:\GitHub\filters\20191029_ephysBoothSpkrRight_300-80k_fs400k.mat';
load(filterName);

%Prepare stimulus sequence
uniqueEvents = length(optoOffset)+2; %+2 for sound w/o laser, and laser w/o sound
index = [1 1 0; 2 0 1; 3 1 1; 4 1 1; 5 1 1; 6 1 1; 7 1 1; 8 1 1; 9 1 1]; %column 2 is sound, 3 is laser
index(:,4) = [0 0 optoOffset]'; %column 4 is laser offset
order = repmat(index(:,1)',[1 repeats]);
order=order(randperm(length(order)));

%Create stimulus
totalSamples = (repeats*uniqueEvents*(burstDuration+IBI)+preStimSilence)*sampleRate/1000;
stimulus = zeros(3,totalSamples);

clickSamples = burstDuration*sampleRate/1000;
IBIsamples = IBI*sampleRate/1000;
preStimSilenceSamples = preStimSilence*sampleRate/1000;
rampSamples = rampDuration*sampleRate/1000;
pulseSamples = optoPulseDuration*sampleRate/1000;

for i = 1:length(order)
    click = rand(1,clickSamples);
    click = applyRamp_AMW(click,rampSamples);
    attn = 70-intensity;
    click = click.*10^(-attn/20);
    
    pulse = zeros(1,pulseSamples*2);
    pulse(1:pulseSamples) = 5;
    pulseCount = ceil(clickSamples/pulseSamples);
    pulse = repmat(pulse,[1 pulseCount]);
    pulse=pulse(1:clickSamples);

    sampleIndex = (i-1)*(clickSamples+IBIsamples)+1+preStimSilenceSamples;
    stimulus(2,sampleIndex:(sampleIndex+clickSamples-1)) = 5;

    if index(order(i),2)==1 %sound present
        stimulus(1,sampleIndex:(sampleIndex+clickSamples-1)) = click;
    end

    if index(order(i),3)==1 %laser present
        offset = index(order(i),4)*sampleRate/1000;
        stimulus(3,(sampleIndex+offset):(sampleIndex+clickSamples-1+offset)) = pulse;
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
stimInfo.burstDuration = burstDuration; % noise duration in ms
stimInfo.optoOffsets = optoOffset;
stimInfo.optoPulseDuration = optoPulseDuration;
stimInfo.index = index;
stimInfo.order = order;
stimInfo.preStimSilence = preStimSilence;
stimInfo.repeats = repeats;
stimInfo.intensity = intensity; % in dB
stimInfo.IBI = IBI; % inter stim interval duration in ms
stimInfo.filter = filterName;

save([fileName '_stimInfo.mat'],'stimInfo');