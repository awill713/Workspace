% Uses the nidaq card for sound instead of PsychPortAudio, therefore
% creates a single long auditory stimulus to relay to the nidaq card
% instead of the sequential buffer loaded to the PsychPort speaker.

sca;
close all;
clear all;

mouse = 'AW000';
recDate = num2str(yyyymmdd(datetime));
saveFolder = fullfile('D:\stimuli\Aaron',mouse,recDate);
filterName = 'D:\GitHub\filters\20191003_ephysBoothSpkr_300-70k_fs400k.mat';
% load(filterName);

%% Input stimulus parameters

% Auditory stimulus parameters
sampleFreq = 40000;
fundFreq = 440;
dbRange = [-2 0];

repeats = 1;
stimDuration = 1; %seconds
isi = 1; %seconds, interstimulus interval

preStimSilence = 5; %seconds
rampDuration = 5; %milliseconds

% Visual stimulus parameters
circleColor = [0 0 0];

circleMinMax = [100 1000]; %pixel size
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
soundSamples = stimDuration*sampleFreq;
rampSamples = rampDuration*sampleFreq/1000;
attn = 0;

loomRamp(1,:) = logspace(dbRange(1),dbRange(2),soundSamples);
loomRamp(2,:) = logspace(dbRange(2),dbRange(1),soundSamples);
loomRamp(3,:) = ones(1,soundSamples);
loomRamp(4,:) = zeros(1,soundSamples);

coef = fundFreq*2*pi/sampleFreq;
soundApproach = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(1,:),rampSamples).*10^(-attn/20);
soundRecede = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(2,:),rampSamples).*10^(-attn/20);
soundStatic = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(3,:),rampSamples).*10^(-attn/20);
soundNone = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(4,:),rampSamples).*10^(-attn/20);

totalStimSeconds = preStimSilence + size(index,1)*repeats*(stimDuration+isi);
totalStimSamples = totalStimSeconds*sampleFreq;
soundStimulus = zeros(2,totalStimSamples);
for event = 1:length(order)
    sampleIndex = eventTimes(event)*sampleFreq + 1;
    if index(order(event),2) == 1
        tempStim = soundApproach;
    elseif index(order(event),2) == 2
        tempStim = soundRecede;
    elseif index(order(event),2) == 3
        tempStim = soundStatic;
    elseif index(order(event),2) == 4
        tempStim = soundNone;
    end
    soundStimulus(1,sampleIndex:(sampleIndex+soundSamples-1)) = tempStim;
    soundStimulus(2,sampleIndex:(sampleIndex+soundSamples-1)) = 5;
end
% soundStimulus(1,:) = conv(soundStimulus(1,:),filterName,'same');

%% Initialize NIDAQ session and load soundStimulus
daqreset;
session = daq.createSession('ni');
session.rate = sampleFreq;
for ch = 1:length(channelsOut)
    addAnalogOutputChannel(session,'dev3',channelsOut(ch),'Voltage');
end
session.queueOutputData(soundStimulus);

%% Setup Psychtoolbox visual stuff

PsychDefaultSetup(2);

screens = Screen('Screens');
screenNumber = max(screens);
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
grey = white / 2;

Screen('Preference', 'SkipSyncTests', 1);
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
circleSizeList(3,:) = ones(1,visStimFrames).*(sum(circleMinMax)/2); %static


%% Present stimulus
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

disp('Presenting stimulus');

session.startBackground();
Priority(topPriorityLevel);

for e = 1:length(lightOn)
    timeOn = lightOn(e);
    timeOff = lightOff(e);
    direction = lightDirectionOrder(e);
    
    for f = 1:visStimFrames
        dotSize = circleSizeList(direction,f);
        circleSize = [0 0 dotSize dotSize];
        centeredRect = CenterRectOnPointd(circleSize,xCenter,yCenter);
        Screen('FillOval', window, circleColor,centeredRect,maxDiameter);
        
        vbl = Screen('Flip',window,triggerTime + timeOn + (f-1)*ifi + 0.5*ifi);

    end
    
    Screen('FillRect',window,grey);
    vbl = Screen('Flip',window,vbl + 0.5*ifi);
    Screen('DrawingFinished',window);
end

disp('Finished stimulus presentation');

sca;

%% Save experiment info
exptInfo.index = index;
exptInfo.order = order;
exptInfo.sampleRate = sampleFreq;
exptInfo.stimDuration = stimDuration;
exptInfo.ISI = isi;
exptInfo.preStimSilence = preStimSilence;
exptInfo.repeats = repeats;
exptInfo.eventTimes = eventTimes;
exptInfo.eventsOut = eventsOut;
exptInfo.fundFreq = fundFreq;
exptInfo.dbRange = dbRange;
exptInfo.circleColor = circleColor;
exptInfo.circleMinMax = circleMinMax;
exptInfo.filter = filterName;
exptInfo.monitorFrameRate = ifi;

fileName = [recDate '_' mouse '_AVloomingComplexTones_exptInfo'];
% save(fullfile(saveFolder,fileName),'exptInfo');
