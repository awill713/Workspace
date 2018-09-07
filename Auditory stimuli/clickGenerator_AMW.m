
clear
seed = 5;
rng(seed)

fileName = ['E:\stimuli\Aaron\WideField\' datestr(now,'yyyymmdd') '_noiseClicks_20rep_10Hz10clicks7sISI100dB_400k_' sprintf('%03d',seed) ]; %don't append .wav or .mat
% filename = ['E:\stimuli\Kath\2P_FRA\' datestr(now,'yyyymmdd') '_pureTones_100ms_5k-32k_70dB_' sprintf('%03d',seed) ];

sampleRate = 400e3;

clickDuration = 25; %ms
clicksPerBurst = 10;
burstDuration = 1000; %ms
burstRepeats = 20;
IBI = 7000; %interburst interval
intensity = 100; %db
rampDuration = 5; %ms

%Load filter
% filterName = 'E:\calibration\Filters\20180427_2Pspkr_IC_NidaqInvFilt_3k-80k_fs400k.mat';
% load(filterName);

%Create stimulus
totalSamples = burstRepeats*(burstDuration+IBI)*sampleRate/1000;
stimulus = zeros(2,totalSamples);

clickSamples = clickDuration*sampleRate/1000;
interClickTime = (burstDuration/clicksPerBurst-clickDuration)/1000; %seconds
interClickSamples = interClickTime*sampleRate;

IBIsamples = IBI*sampleRate/1000;

rampSamples = rampDuration*sampleRate/1000;

for i = 1:burstRepeats
    click = rand(1,clickSamples);
    click = applyRamp_AMW(click,rampSamples);
    attn = 70-intensity;
    click = click.*10^(-attn/20);
    click = [click zeros(1,interClickSamples)];
    
    burst = repmat(click,1,clicksPerBurst);
    
    sampleIndex = (i-1)*(length(burst)+IBIsamples)+1;
    stimulus(1,sampleIndex:(sampleIndex+length(burst)-1)) = burst;
    stimulus(2,sampleIndex:(sampleIndex+length(burst)-1)) = 5;
end
    
% stimulus = conv(stimulus(1,:),'FILT','same');

stimulus = stimulus';

minutes = floor(length(stimulus)/sampleRate/60);
sec = (length(stimulus)/sampleRate/60 - minutes)*60;
disp(['Stimulus is ' num2str(minutes) ' minutes and ' num2str(sec) ' seconds']);

%Save files
audiowrite([fileName '.wav'],(stimulus/10),sampleRate);


stimInfo.seed = seed;
stimInfo.fileName = fileName;
stimInfo.sampleRate = sampleRate;
stimInfo.clickDuration = clickDuration; % noise duration in ms
stimInfo.clicksPerBurst = clicksPerBurst; % number of noise clicks
stimInfo.burstDuration = burstDuration; % in ms
stimInfo.intensity = intensity; % in dB
stimInfo.ISI = IBI; % inter stim interval duration in ms
stimInfo.repeats = burstRepeats;
stimInfo.filter = filterName;

save([fileName '_stimInfo.mat'],'stimInfo')

