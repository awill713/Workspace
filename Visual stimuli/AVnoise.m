
sca;
close all;
clear all;

%% Noise click parameters

sampleFreq = 48000;

repeats = 3;
clickDuration = 1; %seconds
interClickInterval = 2; %seconds

preStimSilence = 5; %seconds

rampDuration = 5; %milliseconds


%% Compile parameters
index = [1 1 0; 2 0 1; 3 1 1];

order = repmat(index(:,1)',[1 repeats]);
order = order(randperm(length(order)));

eventTimes = zeros(1,length(order));
startToStart = interClickInterval + clickDuration;
eventTimes(1) = preStimSilence;
for i = 2:length(order)
    eventTimes(i) = eventTimes(i-1) + startToStart;
end

%% Make noise click, silence, and prestimulus trigger
clickSamples = clickDuration*sampleFreq;
rampSamples = rampDuration*sampleFreq/1000;
click = rand(1,clickSamples);
click = applyRamp_AMW(click,rampSamples);
attn = 0;
click = click.*10^(-attn/20);
click = [click;click];

silence = zeros(2,clickSamples);

preStimSamples = preStimSilence * sampleFreq;
preStimTrigger = zeros(2,preStimSamples);
preStimTrigger(2,1) = 5;

%% Initialize PsychSound and design buffer
InitializePsychSound(1);

waitForDeviceStart = 1;
pahandle = PsychPortAudio('Open',[],3,1,sampleFreq,[2 1]);
PsychPortAudio('Volume',pahandle,1);

totalStimSeconds = preStimSilence + repeats*(clickDuration+interClickInterval);
PsychPortAudio('GetAudioData',pahandle,totalStimSeconds*1.5);

buffer = [];
buffer(1) = PsychPortAudio('CreateBuffer', pahandle, silence);
buffer(2) = PsychPortAudio('CreateBuffer', pahandle, click);
buffer(3) = PsychPortAudio('CreateBuffer', pahandle, preStimTrigger);

PsychPortAudio('UseSchedule', pahandle, 1);
cmdCode = -(1+8);
PsychPortAudio('AddToSchedule',pahandle,buffer(3));
PsychPortAudio('AddToSchedule',pahandle,cmdCode,preStimSilence);

for event = 1:length(order)
    if index(order(event),2) == 1
        PsychPortAudio('AddToSchedule',pahandle,buffer(2));
    elseif index(order(event),2) == 0
        PsychPortAudio('AddToSchedule',pahandle,buffer(1));
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
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,black);
Screen('BlendFunction',window,'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');

ifi = Screen('GetFlipInterval',window);
topPriorityLevel = MaxPriority(window);

clickLengthFrames = round(clickDuration / ifi);
interClickIntervalFrames = round(interClickInterval / ifi);

waitframes = 1;

eventsWithLight = find(index(order,3));
lightOn = eventTimes(eventsWithLight);
lightOff = lightOn + clickDuration;


%% Present stimulus
timeOffset = 2;
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

disp('Presenting stimulus');

startTime = PsychPortAudio('Start', pahandle, [], triggerTime+timeOffset,1);
Priority(topPriorityLevel);
for e = 1:length(lightOn)
    timeOn = lightOn(e);
    timeOff = lightOff(e);
    
    Screen('FillRect',window,white);
    vbl = Screen('Flip',window,startTime + timeOn - 0.5*ifi);
    Screen('DrawingFinished',window);
    
    Screen('FillRect',window,black);
    vbl = Screen('Flip',window,startTime + timeOff - 0.5*ifi);
    Screen('DrawingFinished',window);
end

PsychPortAudio('Stop',pahandle,1,1);
disp('Finished stimulus presentation');

framesOut = PsychPortAudio('GetAudioData',pahandle);

PsychPortAudio('Close',pahandle);
sca;
