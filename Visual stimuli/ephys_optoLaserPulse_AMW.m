
% This script makes a stimulus of optogenetic laser presented in 1 ms and 5
% ms pulses to identify direct projections from the stimulated region to
% the recording region.

clear
seed = 5;
rng(seed)

fileName = ['D:\stimuli\Aaron\Opto\' datestr(now,'yyyymmdd') '_optoLaserPulse_1and5ms_400k_' sprintf('%03d',seed) ]; %don't append .wav or .mat

sampleRate = 400e3;

optoPulseDuration = [1 5]; %ms
repeats = 50;
IBI = 500; %interburst interval, ms
preStimSilence = 1000; %ms

%Prepare stimulus sequence
uniqueEvents = length(optoPulseDuration); %laser on and laser off
index = zeros(uniqueEvents,2); %index,laser duration
for i = 1:size(index,1)
    index(i,1) = i; %index
    index(i,2) = optoPulseDuration(i);
end
order = repmat(index(:,1)',[1 repeats]);
order=order(randperm(length(order)));

%Create stimulus
eventDuration = max(optoPulseDuration); %all events are going to align to this, for ease of event timing
totalSamples = (repeats*uniqueEvents*(eventDuration+IBI)+preStimSilence)*sampleRate/1000;
stimulus = zeros(3,totalSamples);

eventSamples = eventDuration*sampleRate/1000;
IBIsamples = IBI*sampleRate/1000;
preStimSilenceSamples = preStimSilence*sampleRate/1000;

for i = 1:length(order)
    duration = index(order(i));
    pulseSamples = duration*sampleRate/1000;
    
    pulse = zeros(1,eventSamples);
    pulse(1:pulseSamples) = 5;
    
    sampleIndex = (i-1)*(eventSamples+IBIsamples)+1+preStimSilenceSamples;
    stimulus(2,sampleIndex:(sampleIndex+eventSamples-1)) = 5;
    stimulus(3,sampleIndex:(sampleIndex+eventSamples-1)) = pulse;
end

stimulus = stimulus'./10;

minutes = floor(length(stimulus)/sampleRate/60);
sec = (length(stimulus)/sampleRate/60 - minutes)*60;
disp(['Stimulus is ' num2str(minutes) ' minutes and ' num2str(sec) ' seconds']);

%Save files
audiowrite([fileName '.wav'],(stimulus),sampleRate);


stimInfo.seed = seed;
stimInfo.fileName = fileName;
stimInfo.sampleRate = sampleRate;
stimInfo.optoPulseDuration = optoPulseDuration;
stimInfo.repeats = repeats;
stimInfo.index = index;
stimInfo.order = order;
stimInfo.preStimSilence = preStimSilence;
stimInfo.IBI = IBI; % inter stim interval duration in ms

save([fileName '_stimInfo.mat'],'stimInfo');