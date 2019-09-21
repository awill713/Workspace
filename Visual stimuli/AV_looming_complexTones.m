
sca;
close all;
clear all;

%% Input stimulus parameters

% Auditory stimulus parameters
sampleFreq = 40000;
fundFreq = 440;

repeats = 1;
stimDuration = 1; %seconds
isi = 1; %seconds, interstimulus interval

preStimSilence = 5; %seconds
rampDuration = 5; %milliseconds

% Visual stimulus parameters
circleColor = [0 0 0];

circleMinMax = [100 1000]; %pixel size
maxDiameter = max(circleMinMax)*1.01;


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

loomRamp(1,:) = logspace(-2,0,soundSamples);
loomRamp(2,:) = logspace(0,-2,soundSamples);
loomRamp(3,:) = ones(1,soundSamples);
loomRamp(4,:) = zeros(1,soundSamples);

coef = fundFreq*2*pi/sampleFreq;
soundApproach = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(1,:),rampSamples).*10^(-attn/20);
soundRecede = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(2,:),rampSamples).*10^(-attn/20);
soundStatic = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(3,:),rampSamples).*10^(-attn/20);
soundNone = applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(4,:),rampSamples).*10^(-attn/20);

preStimSamples = preStimSilence * sampleFreq;
preStimTrigger = zeros(1,preStimSamples);


%% Initialize PsychSound and design buffer
InitializePsychSound(1);

waitForDeviceStart = 1;
pahandle = PsychPortAudio('Open',2,1,1,sampleFreq,1);
PsychPortAudio('Volume',pahandle,1);

buffer = [];
buffer(1) = PsychPortAudio('CreateBuffer', pahandle, soundApproach);
buffer(2) = PsychPortAudio('CreateBuffer', pahandle, soundRecede);
buffer(3) = PsychPortAudio('CreateBuffer', pahandle, soundStatic);
buffer(4) = PsychPortAudio('CreateBuffer', pahandle, soundNone);
buffer(5) = PsychPortAudio('CreateBuffer', pahandle, preStimTrigger);

PsychPortAudio('UseSchedule', pahandle, 1);
cmdCode = -(1+8);
PsychPortAudio('AddToSchedule',pahandle,buffer(5));
PsychPortAudio('AddToSchedule',pahandle,cmdCode,preStimSilence);

for event = 1:length(order)
    if index(order(event),2) == 1
        PsychPortAudio('AddToSchedule',pahandle,buffer(1));
    elseif index(order(event),2) == 2
        PsychPortAudio('AddToSchedule',pahandle,buffer(2));
    elseif index(order(event),2) == 3
        PsychPortAudio('AddToSchedule',pahandle,buffer(3));
    elseif index(order(event),2) == 4
        PsychPortAudio('AddToSchedule',pahandle,buffer(4));
    end
    PsychPortAudio('AddToSchedule',pahandle,cmdCode,startToStart);
end

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
circleSizeList = linspace(circleMinMax(1),circleMinMax(2),visStimFrames);


%% Present stimulus
timeOffset = 2;
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

disp('Presenting stimulus');

startTime = PsychPortAudio('Start', pahandle, [], triggerTime+timeOffset,1);
Priority(topPriorityLevel);

for e = 1:length(lightOn)
    eventsWithLight(e)
    timeOn = lightOn(e);
    timeOff = lightOff(e);
    direction = lightDirectionOrder(e);
    
    time = 0;
    
    for f = 1:visStimFrames
        dotSize = circleSizeList(f);
        circleSize = [0 0 dotSize dotSize];
        centeredRect = CenterRectOnPointd(circleSize,xCenter,yCenter);
        Screen('FillOval', window, circleColor,centeredRect,maxDiameter);
        
        vbl = Screen('Flip',window,startTime + timeOn + (f-1)*ifi + 0.5*ifi);
        
%         time sca= time+ifi;
    end
    
    Screen('FillRect',window,grey);
    vbl = Screen('Flip',window,vbl + 0.5*ifi);
    Screen('DrawingFinished',window);
end

PsychPortAudio('Stop',pahandle,1,1);
disp('Finished stimulus presentation');

PsychPortAudio('Close',pahandle);
sca;

%% Setup Psychtoolbox visual stuff
% PsychDefaultSetup(2);
% 
% screens = Screen('Screens');
% screenNumber = max(screens);
% black = BlackIndex(screenNumber);
% white = WhiteIndex(screenNumber);
% grey = white / 2;
% 
% Screen('Preference', 'SkipSyncTests', 1);
% [window, windowRect] = PsychImaging('OpenWindow',screenNumber,grey);
% Screen('BlendFunction',window,'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');
% 
% [screenXpixels, screenYpixels] = Screen('WindowSize',window);
% [xCenter, yCenter] = RectCenter(windowRect);
% 
% ifi = Screen('GetFlipInterval',window);
% 
% vbl = Screen('Flip',window);
% waitframes = 1;
% 
% topPriorityLevel = MaxPriority(window);
% Priority(topPriorityLevel);
% 
% while ~KbCheck
%     
%     dotSize = round(amplitude * sin(frequency * time * 2*pi) + 1) + amplitude + sizeOffset;
%     circleSize = [0 0 dotSize dotSize];
%     centeredRect = CenterRectOnPointd(circleSize,xCenter,yCenter);
%     Screen('FillOval', window, circleColor,centeredRect,maxDiameter);
%     
%     vbl = Screen('Flip',window,vbl+(waitframes - 0.5)*ifi);
%     
%     time = time+ifi;
% end
% 
% sca;

