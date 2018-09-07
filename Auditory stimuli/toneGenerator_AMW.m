
clear
seed = 5;
rng(seed)

fileName = ['E:\stimuli\Aaron\2P_FRA\' datestr(now,'yyyymmdd') '_attnTones_50-70db_100ms_4k-70k_' sprintf('%03d',seed) ]; %don't append .wav or .mat
% fileName = ['E:\stimuli\Aaron\2P_FRA\' datestr(now,'yyyymmdd') '_pureTones_100ms_4k-70k_70dB_' sprintf('%03d',seed) ]; %don't append .wav or .mat
% filename = ['E:\stimuli\Kath\2P_FRA\' datestr(now,'yyyymmdd') '_pureTones_100ms_5k-32k_70dB_' sprintf('%03d',seed) ];

sampleRate = 400e3;

%Design stimulus parameters
freqLow = 440;
freqHigh = 880;
uniqueTones = 5;

repeats = 5;

attenuations = [0];

toneDuration = 1000; %milliseconds
ITI = 3000; %milliseconds
rampDuration = 5; %milliseconds

%Load filter
% filterName = 'E:\calibration\Filters\20180427_2Pspkr_IC_NidaqInvFilt_3k-80k_fs400k.mat';
% load(filterName);

%Create stimulus sequence
logToneRange = linspace(log10(freqLow),log10(freqHigh),uniqueTones);

index = zeros(uniqueTones*length(attenuations),3);
for i = 1:size(index,1)
    index(i,1) = i;
    index(i,2) = 10^logToneRange(floor((i-1)/3)+1);
    index(i,3) = attenuations(mod(i-1,1)+1);
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

stimulus = zeros(2,length(toneOrder)*presentationSamples);

for i = 1:length(toneOrder)
        tempFreq = index(toneOrder(i),2);
        coef = tempFreq*2*pi/sampleRate;
        
        tempSignal(1:toneSamples) = cos(coef*(1:toneSamples));
        
        tempSignal = applyRamp_AMW(tempSignal,rampSamples);
        
        tempAtten = 10^(index(toneOrder(i),3)/20);
        tempSignal = tempSignal.*tempAtten;
        
        sampleIndex = (i-1)*presentationSamples+1;
        stimulus(1,sampleIndex:(sampleIndex+toneSamples-1)) = tempSignal;
        stimulus(2,sampleIndex:(sampleIndex+toneSamples-1)) = 5;
end

% stimulus(1,:) = conv(stimulus(1,:),'FILT','same');

stimulus = stimulus';

minutes = floor(length(stimulus)/sampleRate/60);
sec = (length(stimulus)/sampleRate/60 - minutes)*60;
disp(['Stimulus is ' num2str(minutes) ' minutes and ' num2str(sec) ' seconds']);

%Save files
% audiowrite([fileName '.wav'],(stimulus/10),sampleRate);
audiowrite('SampleStimulus.wav',(stimulus),sampleRate);

% stimInfo.fileName = fileName;
% stimInfo.sampleRate = sampleRate;
% stimInfo.frequencies = 10.^logToneRange;
% stimInfo.uniqueTones = uniqueTones;
% stimInfo.attenuations = attenuations;
% stimInfo.toneDuration = toneDuration; % tone duration in ms
% stimInfo.ITI = ITI; % inter tone interval duration in ms
% stimInfo.repeats = repeats; % number of repeats of each tone
% stimInfo.filterName = filterName;
% stimInfo.order = toneOrder;
% stimInfo.index = index;
% stimInfo.stimGenFunc = 'toneGenerator.m';
% 
% save([fileName '_stimInfo.mat'],'stimInfo')
