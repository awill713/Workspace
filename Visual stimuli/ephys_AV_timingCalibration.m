
close all;
clear all;

daqreset;

%% Design stuff
stimDur = 1; %seconds
isi = 1; %seconds
preStimSilence = 5; %seconds
repeats = 100;

sampleRate = 400e3;
rampDuration = 5; %milliseconds

soundInputChannel = 1;
lightInputChannel = [];
soundOutputChannel = 0;

filterPath = 'D:\GitHub\filters\20191008_ephysBoothSpkr_300-70k_fs400k_TESTACTUAL';
% load(filterPath);

%% Generate auditory stimulus
stimulusSamples = repeats*(stimDur + isi)*sampleRate + preStimSilence;
soundStimulus = zeros(1,stimulusSamples);
eventTimeSeconds = [0:repeats-1]*(stimDur+isi) + preStimSilence;
eventTimeSamples = eventTimeSeconds*sampleRate;
rampSamples = rampDuration*sampleRate/1000;

noiseSamples = stimDur*sampleRate;
noise = applyRamp_AMW(randn(1,noiseSamples),rampSamples);
for r = 1:repeats
    sampleStart = eventTimeSamples(r)+1;
    soundStimulus(sampleStart:sampleStart+noiseSamples-1) = noise;
end

% soundStimulus = conv(soundStimulus,FILT,'same');

%% Set up nidaq stuff
session = daq.createSession('ni');
session.Rate = sampleRate;
addAnalogOutputChannel(session,'dev3',soundOutputChannel,'Voltage');
addAnalogInputChannel(session,'dev3',[soundInputChannel lightInputChannel],'Voltage');
session.Channels(2).InputType = 'SingleEnded';

session.queueOutputData(soundStimulus');

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

lightOn = eventTimeSeconds;
lightOff = lightOn + stimDur;

visStimFrames = round(stimDuration/ifi);
circleSizeList(1,:) = linspace(circleMinMax(1),circleMinMax(2),visStimFrames); %approaching
circleSizeList(2,:) = linspace(circleMinMax(2),circleMinMax(1),visStimFrames); %receding
circleSizeList(3,:) = ones(1,visStimFrames).*(sum(circleMinMax)/2); %static

%% Present stimulus
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

disp('Presenting stimulus');

[soundInput lightInput] = session.startBackground();
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

disp('Finished stimulus presentation');

sca;

%% Analyze timing

