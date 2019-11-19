
sca;
close all;
clear all;
daqreset;

mouse = 'AW102';
recDate = num2str(yyyymmdd(datetime));
saveFolder = fullfile('D:\data\',mouse,recDate);
filterName = 'D:\GitHub\filters\20191029_ephysBoothSpkrRight_300-80k_fs400k';
load(filterName);

%% Input stimulus parameters

% Auditory stimulus parameters
sampleRate = 400e3;

avOffset = 0; %milliseconds

repeats = 20;
stimDuration = 0.1; %seconds
isi = 0.5; %seconds, interstimulus interval

intensities = [0 40 50 60 70]; %db

preStimSilence = 5; %seconds
rampDuration = 5; %milliseconds

%Nidaq settings
channelsOut = [0 1]; %speaker and events, respectively
soundInputChannel = 1;
lightInputChannel = 5;
soundOutputChannel = 0;

%% Generate stimulus indexing, order, and event times
uniqueEvents = length(intensities)*2; %light on and light off
index = zeros(uniqueEvents,3); %index,intensity,light
for i = 1:size(index,1)
    index(i,1) = i; %index
    index(i,2) = intensities(mod(i-1,length(intensities))+1); %intensity
    index(i,3) = floor((i-1)/length(intensities)); %light
end
order = repmat(index(:,1)',[1 repeats]);
order=order(randperm(length(order)));

eventTimes = zeros(1,length(order)); %eventTimes is unused but can give a rough estimate if for some reason event TTLs are unregistered/unreliable
startToStart = isi + stimDuration;
eventTimes(1) = preStimSilence;
for i = 2:length(order)
    eventTimes(i) = eventTimes(i-1) + startToStart;
end

%% Make sound stimuli
soundSamples = stimDuration*sampleRate;
rampSamples = rampDuration*sampleRate/1000;

isiVector = zeros(1,isi*sampleRate-avOffset*sampleRate/1000);
eventVector = [zeros(1,avOffset*sampleRate/1000) ones(1,stimDuration*sampleRate)*5 isiVector];

soundDatabase{1} = [zeros(1,(stimDuration+isi)*sampleRate); eventVector];
for i = 2:length(intensities)
    intensity = intensities(i);
    noise = randn(1,soundSamples);
    noise = 10^((intensity-70)/20) * noise;
    noise = [zeros(1,avOffset*sampleRate/1000) applyRamp_AMW(noise,rampSamples) isiVector];
%     noise = [conv(noise,FILT,'same'); eventVector];
    noise = [noise; eventVector];
    soundDatabase{i} = noise;
end

%% Make TTL to tell Cheetah the stimulus is beginning, arbitrary length

ttlLength = 500; %ms
initialTTL = [zeros(1,ttlLength/2*sampleRate/1000) ones(1,ttlLength*sampleRate/1000)*5 zeros(1,ttlLength/2*sampleRate/1000)];
initialTTL = [zeros(1,length(initialTTL)); initialTTL];

%% Initialize NIDAQ session and load soundStimulus
session = daq.createSession('ni');
session.Rate = sampleRate;

addAnalogOutputChannel(session,'dev3',channelsOut,'Voltage');
addAnalogInputChannel(session,'dev3',[soundInputChannel lightInputChannel],'Voltage');
session.Channels(2).InputType = 'SingleEnded';
session.Channels(3).InputType = 'SingleEnded';

fid1 = fopen('log.bin','w');
session.addlistener('DataAvailable',@(src,event) logData(src,event,fid1));

%% Setup Psychtoolbox visual stuff

PsychDefaultSetup(2);

screens = Screen('Screens');
screenNumber = max(screens);
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
grey = white / 2;

% Screen('Preference', 'SkipSyncTests', 1);
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,black);
Screen('BlendFunction',window,'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');

[screenXpixels, screenYpixels] = Screen('WindowSize',window);

ifi = Screen('GetFlipInterval',window);
topPriorityLevel = MaxPriority(window);


%% Present stimulus
session.queueOutputData(initialTTL');
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

session.startBackground();
disp('Presenting stimulus');

pause(preStimSilence);

for evNum = 1:repeats
    type = 10;
    soundFile = soundDatabase{intensities==index(type,2)};

    lightPresent = index(type,3);
    
    session.queueOutputData(soundFile');
    wait(session);
    
    if lightPresent
        session.startBackground();
        
        Screen('FillRect',window,white);
        vbl = Screen('Flip',window);
        Screen('DrawingFinished',window);
    
        Screen('FillRect',window,black);
        vbl = Screen('Flip',window,vbl + stimDuration - 0.5*ifi);
        Screen('DrawingFinished',window);
    else
        session.startBackground();
    end
    
    wait(session);
end
    
disp('Finished stimulus presentation');
sca;


%% Analyze timing

fid2 = fopen('log.bin','r');
[data count] = fread(fid2,[3,inf],'double');
fclose(fid2);

soundInput = data(2,:);
[fb, fa] = butter(5, 2*300 / sampleRate, 'high');
soundInput = filter(fb,fa,data(2,:));
lightInput = data(3,:);

figure; hold on;
plot(soundInput);
plot(lightInput);

baselineSoundVolt = max(soundInput((end-sampleRate*isi*0.2):end));
baselineSoundVolt = max(soundInput(1:0.3*isi*sampleRate));
soundVoltThresh = baselineSoundVolt*50;
soundTimeThresh = 0.5*sampleRate*isi;
thresh = soundInput>soundVoltThresh;
supraThresh = find(thresh==1);
diffSupra = diff(supraThresh);
longDiff = find(diffSupra>soundTimeThresh); 
soundOnset = supraThresh([1 longDiff+1]);

baselineLightVolt = max(soundInput(1:1.5e6));
% lightVoltThresh = baselineLightVolt*1.5;
lightVoltThresh = 2.1;
lightTimeThresh = 0.2*sampleRate*isi;
thresh = lightInput>lightVoltThresh;
supraThresh = find(thresh==1);
diffSupra = diff(supraThresh);
longDiff = find(diffSupra>lightTimeThresh);
lightOnset = supraThresh([1 longDiff+1]);


offset = (soundOnset - lightOnset)./sampleRate *1000; %milliseconds
figure;histogram(offset,'BinWidth',1)
figure;plot(offset)
mean(offset)


