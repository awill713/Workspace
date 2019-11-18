% Uses the nidaq card for sound instead of PsychPortAudio. Timing between
% the two is pretty bad, so can't use a single long auditory stimulus.
% Instead, present individual sound/light trials, telling Psychtoolbox to
% immediately start flipping screen after sound stimulus starts
% (session.startBackground()). Incorporate a delay in the sound stimulus in
% order to account for the delay in starting the visual stimulus relative
% to sound stimulus using this method, usually ~40-50ms. This stimulus
% seems to be ~41ms.


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
dbRange = [50 70];

avOffset = 41; %milliseconds

repeats = 20;
stimDuration = 1; %seconds
isi = 3; %seconds, interstimulus interval

preStimSilence = 5; %seconds
rampDuration = 5; %milliseconds

% Visual stimulus parameters
circleColor = [0 0 0];

circleMinMax = [50 700]; %pixel size
maxDiameter = max(circleMinMax)*1.01;

%Nidaq settings
channelsOut = [0 1]; %speaker and events, respectively

%% Generate stimulus indexing, order, and event times
index(:,1) = 1:16;
index(:,2) = ceil([1:16]./4)'; %1=approach, 2=recede,3=static,4=nothing
index(:,3) = repmat([1 2 3 4],[1 4])';

order = repmat(index(:,1)',[1 repeats]);
order = order(randperm(length(order)));

eventTimes = zeros(1,length(order));
startToStart = isi + stimDuration;
eventTimes(1) = preStimSilence;
for i = 2:length(order)
    eventTimes(i) = eventTimes(i-1) + startToStart;
end

%% Make sound stim (approaching, receding, silence, prestimulus trigger)
soundSamples = stimDuration*sampleRate;
rampSamples = rampDuration*sampleRate/1000;

loomRamp(1,:) = logspace((dbRange(1)-70)/20,(dbRange(2)-70)/20,soundSamples);
loomRamp(2,:) = logspace((dbRange(2)-70)/20,(dbRange(1)-70)/20,soundSamples);
loomRamp(3,:) = ones(1,soundSamples);
loomRamp(4,:) = zeros(1,soundSamples);

isiVector = zeros(1,isi*sampleRate-avOffset*sampleRate/1000);
eventVector = [zeros(1,avOffset*sampleRate/1000) ones(1,stimDuration*sampleRate)*5 isiVector];

soundApproach = [zeros(1,avOffset*sampleRate/1000) applyRamp_AMW(randn(1,soundSamples).*loomRamp(1,:),rampSamples) isiVector];
soundRecede = [zeros(1,avOffset*sampleRate/1000) applyRamp_AMW(randn(1,soundSamples).*loomRamp(2,:),rampSamples) isiVector];
soundStatic = [zeros(1,avOffset*sampleRate/1000) applyRamp_AMW(randn(1,soundSamples).*loomRamp(3,:),rampSamples) isiVector];
soundNone = [zeros(1,avOffset*sampleRate/1000) applyRamp_AMW(randn(1,soundSamples).*loomRamp(4,:),rampSamples) isiVector];

% soundApproach = [conv(soundApproach,FILT,'same'); eventVector];
% soundRecede = [conv(soundRecede,FILT,'same'); eventVector];
% soundStatic = [conv(soundStatic,FILT,'same'); eventVector];
% soundNone = [conv(soundNone,FILT,'same'); eventVector];

soundApproach = [soundApproach; eventVector];
soundRecede = [soundRecede; eventVector];
soundStatic = [soundStatic; eventVector];
soundNone = [soundNone; eventVector];

soundDatabase = {soundApproach;soundRecede;soundStatic;soundNone};

%% Make TTL to tell Cheetah the stimulus is beginning, arbitrary length

ttlLength = 500; %ms
initialTTL = [zeros(1,ttlLength/2*sampleRate/1000) ones(1,ttlLength*sampleRate/1000)*5 zeros(1,ttlLength/2*sampleRate/1000)];
initialTTL = [zeros(1,length(initialTTL)); initialTTL];

%% Initialize NIDAQ session and load soundStimulus
session = daq.createSession('ni');
session.Rate = sampleRate;

addAnalogOutputChannel(session,'dev3',channelsOut,'Voltage');

%% Setup Psychtoolbox visual stuff

PsychDefaultSetup(2);

screens = Screen('Screens');
screenNumber = max(screens);
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
grey = white / 2;

% Screen('Preference', 'SkipSyncTests', 1);
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,grey);
Screen('BlendFunction',window,'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');

[screenXpixels, screenYpixels] = Screen('WindowSize',window);
[xCenter, yCenter] = RectCenter(windowRect);

ifi = Screen('GetFlipInterval',window);
topPriorityLevel = MaxPriority(window);

waitframes = 1;

eventsWithLight = find(index(order,3)~=4);
lightOn = eventTimes(eventsWithLight);
lightOff = lightOn + stimDuration;
lightDirectionOrder = index(order(eventsWithLight)',3);

visStimFrames = round(stimDuration/ifi);
circleSizeList(1,:) = linspace(circleMinMax(1),circleMinMax(2),visStimFrames); %approaching
circleSizeList(2,:) = linspace(circleMinMax(2),circleMinMax(1),visStimFrames); %receding
circleSizeList(3,:) = ones(1,visStimFrames).*(mean(circleMinMax)); %static
circleSizeList(4,:) = ones(1,visStimFrames).*(mean(circleMinMax)); %static, but will be same color as background
circleColorList = {circleColor;circleColor;circleColor;grey};


%% Present stimulus
session.queueOutputData(initialTTL');
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

session.startBackground();
disp('Presenting stimulus');

pause(preStimSilence);

for evNum = 1:length(order)
    type = order(evNum);
    soundFile = soundDatabase{index(type,2)};
    
    color = circleColorList{index(type,3)};
    
    session.queueOutputData(soundFile');
    wait(session);
    
    session.startBackground();
    
    for f = 1:visStimFrames
        dotSize = circleSizeList(index(type,3),f);
        circleSize = [0 0 dotSize dotSize];
        centeredRect = CenterRectOnPointd(circleSize,xCenter,yCenter);
        Screen('FillOval', window, color,centeredRect,maxDiameter);
        
        vbl = Screen('Flip',window);
    end
    
    Screen('FillRect',window,grey);
    vbl = Screen('Flip',window);
    Screen('DrawingFinished',window);
    
    wait(session);
end
    
disp('Finished stimulus presentation');
sca;


%% Save experiment info
stimInfo.index = index;
stimInfo.order = order;
stimInfo.sampleRate = sampleRate;
stimInfo.stimDuration = stimDuration;
stimInfo.ISI = isi;
stimInfo.preStimSilence = preStimSilence;
stimInfo.repeats = repeats;
stimInfo.eventTimes = eventTimes;
stimInfo.dbRange = dbRange;
stimInfo.circleColor = circleColor;
stimInfo.circleMinMax = circleMinMax;
stimInfo.filter = filterName;
stimInfo.monitorFrameRate = ifi;

fileName = [recDate '_' mouse '_AVloomingWhiteNoise_stimInfo'];

if ~exist(saveFolder, 'dir')
   mkdir(saveFolder)
end
save(fullfile(saveFolder,fileName),'stimInfo');