
close all;
clear all;

daqreset;

%% Design stuff
stimDur = 4; %seconds
isi = 1; %seconds
preStimSilence = 5; %seconds
repeats = 100;

sampleRate = 400e3;
rampDuration = 5; %milliseconds

avOffset = 0; %milliseconds

soundInputChannel = 1;
lightInputChannel = 5;
soundOutputChannel = 0;

filterPath = 'D:\GitHub\filters\20191008_ephysBoothSpkr_300-70k_fs400k_TESTACTUAL';
% load(filterPath);

%% Generate auditory stimulus
eventTimeSeconds = [0:repeats-1]*(stimDur+isi) + preStimSilence;
eventTimeSamples = eventTimeSeconds*sampleRate;
rampSamples = rampDuration*sampleRate/1000;

noiseSamples = stimDur*sampleRate;
noise = applyRamp_AMW(randn(1,noiseSamples),rampSamples);

soundStimulus = [zeros(1,avOffset*sampleRate/1000) noise zeros(1,isi*sampleRate-avOffset*sampleRate/1000)];
% soundStimulus = [noise zeros(1,isi*sampleRate)];

% soundStimulus = conv(soundStimulus,FILT,'same');

%% Set up nidaq stuff, saving incoming data to fid1
session = daq.createSession('ni');
session.Rate = sampleRate;
addAnalogOutputChannel(session,'dev3',soundOutputChannel,'Voltage');
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
[xCenter, yCenter] = RectCenter(windowRect);

ifi = Screen('GetFlipInterval',window);
topPriorityLevel = MaxPriority(window);

lightOn = eventTimeSeconds;
lightOff = lightOn + stimDur;
% 
visStimFrames = round(stimDur/ifi);
% circleSizeList(1,:) = linspace(circleMinMax(1),circleMinMax(2),visStimFrames); %approaching
% circleSizeList(2,:) = linspace(circleMinMax(2),circleMinMax(1),visStimFrames); %receding
% circleSizeList(3,:) = ones(1,visStimFrames).*(sum(circleMinMax)/2); %static

%% Present stimulus
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

for r = 1:repeats
    session.queueOutputData(soundStimulus');
    wait(session);
    
    session.startBackground();
    
    Screen('FillRect',window,white);
    vbl = Screen('Flip',window);
    Screen('DrawingFinished',window);
    
    Screen('FillRect',window,black);
    vbl = Screen('Flip',window,vbl+stimDur);
    Screen('DrawingFinished',window);
    
    wait(session);
    
end

disp('Finished stimulus presentation');
fclose(fid1);
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


