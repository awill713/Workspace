% Uses the nidaq card for sound instead of PsychPortAudio. Timing between
% the two is pretty bad, so can't use a single long auditory stimulus.
% Instead, present individual sound/light trials, telling Psychtoolbox to
% immediately start flipping screen after sound stimulus starts
% (session.startBackground()). Incorporate a delay in the sound stimulus in
% order to account for the delay in starting the visual stimulus relative
% to sound stimulus using this method, usually ~40-50ms. This stimulus
% seems to be ~41ms

sca;
close all;
clear all;
daqreset;

mouse = 'AW000';
recDate = num2str(yyyymmdd(datetime));
saveFolder = fullfile('D:\stimuli\Aaron',mouse,recDate);
filterName = 'D:\GitHub\filters\20191008_ephysBoothSpkr_300-70k_fs400k_TESTACTUAL';
load(filterName);

%% Input stimulus parameters

% Auditory stimulus parameters
sampleFreq = 400e3;
fundFreq = 440;
dbRange = [60 80];

repeats = 3;
stimDuration = 1; %seconds
isi = 3; %seconds, interstimulus interval

preStimSilence = 5; %seconds
rampDuration = 5; %milliseconds

% Visual stimulus parameters
circleColor = [0 0 0];

circleMinMax = [100 1000]; %pixel size
maxDiameter = max(circleMinMax)*1.01;

%Nidaq settings
channelsOut = [0 1]; %speaker and events, respectively
soundInputChannel = 1;
lightInputChannel = 5;

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

loomRamp(1,:) = logspace((dbRange(1)-70)/20,(dbRange(2)-70)/20,soundSamples);
loomRamp(2,:) = logspace((dbRange(2)-70)/20,(dbRange(1)-70)/20,soundSamples);
loomRamp(3,:) = ones(1,soundSamples);
loomRamp(4,:) = zeros(1,soundSamples);

isiVector = zeros(1,isi*sampleFreq);
eventVector = [zeros(1,stimDuration*sampleFreq)*5 isiVector];

coef = fundFreq*2*pi/sampleFreq;
soundApproach = [applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(1,:),rampSamples).*10^(-attn/20) isiVector];
soundRecede = [applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(2,:),rampSamples).*10^(-attn/20) isiVector];
soundStatic = [applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(3,:),rampSamples).*10^(-attn/20) isiVector];
soundNone = [applyRamp_AMW(sawtooth(coef*(1:soundSamples),0.5).*loomRamp(4,:),rampSamples).*10^(-attn/20) isiVector];

% soundApproach = [conv(soundApproach,FILT,'same')/10; eventVector];
% soundRecede = [conv(soundRecede,FILT,'same')/10; eventVector];
% soundStatic = [conv(soundStatic,FILT,'same')/10; eventVector];
% soundNone = [conv(soundNone,FILT,'same')/10; eventVector];

soundApproach = [soundApproach; eventVector];
soundRecede = [soundRecede; eventVector];
soundStatic = [soundStatic; eventVector];
soundNone = [soundNone; eventVector];

soundDatabase = {soundApproach;soundRecede;soundStatic;soundNone};

%% Initialize NIDAQ session and load soundStimulus
session = daq.createSession('ni');
session.Rate = sampleFreq;

addAnalogOutputChannel(session,'dev3',channelsOut,'Voltage');
addAnalogInputChannel(session,'dev3',[soundInputChannel lightInputChannel],'Voltage');
session.Channels(3).InputType = 'SingleEnded';
session.Channels(4).InputType = 'SingleEnded';

fid1 = fopen('log.bin','w');
session.addlistener('DataAvailable',@(src,event) logData(src,event,fid1));


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
circleSizeList(3,:) = ones(1,visStimFrames).*(mean(circleMinMax)); %static
circleSizeList(4,:) = ones(1,visStimFrames).*(mean(circleMinMax)); %static, but will be same color as background
circleColorList = {circleColor;circleColor;circleColor;grey};


%% Present stimulus
disp('Hit a key to begin stimulus presentation');
triggerTime = KbStrokeWait;

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
fclose(fid1);
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
exptInfo.fundFreq = fundFreq;
exptInfo.dbRange = dbRange;
exptInfo.circleColor = circleColor;
exptInfo.circleMinMax = circleMinMax;
exptInfo.filter = filterName;
exptInfo.monitorFrameRate = ifi;

fileName = [recDate '_' mouse '_AVloomingComplexTones_exptInfo'];
% save(fullfile(saveFolder,fileName),'exptInfo');

%% Analyze timing

fid2 = fopen('log.bin','r');
[data count] = fread(fid2,[3,inf],'double');
fclose(fid2);

soundInput = data(2,:);
[fb, fa] = butter(5, 2*300 / sampleFreq, 'high');
soundInput = filter(fb,fa,data(2,:));
lightInput = data(3,:);

figure; hold on;
plot(soundInput);
plot(lightInput);


baselineSoundVolt = max(soundInput((end-sampleFreq*isi*00.5):end));
soundVoltThresh = baselineSoundVolt*10;
soundTimeThresh = 0.5*sampleFreq*isi;
thresh = soundInput>soundVoltThresh;
supraThresh = find(thresh==1);
diffSupra = diff(supraThresh);
longDiff = find(diffSupra>soundTimeThresh); 
soundOnset = supraThresh([1 longDiff+1]);

baselineLightVolt = max(soundInput(1:1.5e6));
% lightVoltThresh = baselineLightVolt*1.5;
lightVoltThresh = 0.7;
lightTimeThresh = 0.2*sampleFreq*isi;
thresh = lightInput<lightVoltThresh;
supraThresh = find(thresh==1);
diffSupra = diff(supraThresh);
longDiff = find(diffSupra>lightTimeThresh);
lightOnset = supraThresh([1 longDiff+1]);

evSound = find(index(order,2)~=4); noSound = find(index(order,2)==4);
evLight = find(index(order,3)~=4); noLight = find(index(order,3)==4);
[~, soundNoLightIndx] = intersect(evSound,noLight);
[~, lightNoSoundIndx] = intersect(evLight,noSound);
soundOnset(soundNoLightIndx') = [];
lightOnset(lightNoSoundIndx') = [];
offset = (soundOnset - lightOnset)./sampleFreq *1000; %milliseconds
figure;histogram(offset,'BinWidth',1)
figure;plot(offset)
mean(offset)

