
sca;
close all;
clear all;

%% Create auditory stimulus and events

seed = 5;
rng(seed)

sampleRate = 48000;

%Design stimulus parameters
freqMin = 220;
freqMax = 880;
uniqueFreq = 3;
repeats = 2;
preStimulusSilence = 1; %seconds

toneDuration = 1000; %milliseconds
ITI = 1000; %milliseconds
rampDuration = 5; %milliseconds

%Load filter
% filterName = 'E:\calibration\Filters\20180427_2Pspkr_IC_NidaqInvFilt_3-80k_fs400k.mat';
% load(filterName);

%Create stimulus sequence
freqRange = 10.^linspace(log10(freqMin),log10(freqMax),uniqueFreq);

toneOrder = [];
for r = 1:repeats
    toneOrder = [toneOrder randperm(uniqueFreq)];
end

index = zeros(uniqueFreq,2);
for i = 1:size(index,1)
    index(i,1) = i;
    index(i,2) = freqRange(i);
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
        
        sampleIndex = (i-1)*presentationSamples+1;
        stimulus(1,sampleIndex:(sampleIndex+toneSamples-1)) = tempSignal;
        stimulus(2,sampleIndex:(sampleIndex+toneSamples-1)) = 5;
end

% stimulus(1,:) = conv(stimulus(1,:),FILT,'same');

preStimSilenceSamples = zeros(2,preStimulusSilence*sampleRate);
preStimSilenceSamples(2,1) = 5;

stimulus = [preStimSilenceSamples stimulus];


% stimulus = stimulus';

totalSeconds = length(stimulus)/sampleRate;
minutes = floor(totalSeconds/60);
sec = (totalSeconds/60 - minutes)*60;
disp(['Stimulus is ' num2str(minutes) ' minutes and ' num2str(sec) ' seconds']);

%Save files
% audiowrite([fileName '.wav'],(stimulus/10),sampleRate);
% 
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

%% Initialize PsychPortAudio etc

InitializePsychSound(1);

nrchannels = 2;

startCue = 0; 


waitForDeviceStart = 1;
pahandle = PsychPortAudio('Open',[],3,1,sampleRate,[2 1]);
PsychPortAudio('Volume',pahandle,1);
PsychPortAudio('FillBuffer',pahandle,[stimulus(1,:);stimulus(1,:)]);


PsychPortAudio('GetAudioData',pahandle,2*totalSeconds);
 
disp('\nHit a key to begin stimulus presentation\n\n');
triggerTime = KbStrokeWait;

a = GetSecs; 
PsychPortAudio('Start',pahandle,1,triggerTime,waitForDeviceStart);
disp('Presenting stimulus\n\n');

PsychPortAudio('Stop',pahandle,1,1);
b = GetSecs;
disp('Finished stimulus presentation\n\n');

framesOut = PsychPortAudio('GetAudioData',pahandle);

PsychPortAudio('Close',pahandle);
sca;

% %Screen setup
% PsychDefaultSetup(2);
% 
% screens = Screen('Screens');
% screenNumber = max(screens);
% 
% black = BlackIndex(screenNumber);
% white = WhiteIndex(screenNumber);
% grey = white/2;
% 
% % Screen('Preference', 'SkipSyncTests', 1);
% [window, windowRect] = PsychImaging('OpenWindow',screenNumber,grey);
% Screen('BlendFunction',window,'GL_SRC_ALPHA','GL_ONE_MINUS_SRC_ALPHA');
% [screenXpixels, screenYpixels] = Screen('WindowSize',window);
% 
% ifi = Screen('GetFlipInterval',window);
% [xcenter,yCenter] = RectCenter(windowRect);
% Screen('TextSize',window,70);
% 
% beepLengthFrames = round(beepLengthSecs / ifi);
% beepPauseLengthFrames = round(beepPauseTime / ifi);
% 
% lum = cos(envelopeFreq*2*pi*ifi*(1:beepLengthFrames))/2 + 0.5;
% waitframes = 1;
% 
% tstart = tic;
% for i = 1:beepPauseLengthFrames
%     time = toc(tstart);
%     DrawFormattedText(window,['SILENCE #1'],'center','center',[0 0 1]);
%     Screen('Flip',window);
% end
% PsychPortAudio('Start',pahandle,repetitions,startCue,waitForDeviceStart);
% 
% vbl = Screen('Flip', window);
% for i = 1:beepLengthFrames
% %     DrawFormattedText(window,'BEEP #1','center','center',[0 0 1]);
% %     Screen('Flip',window);
%     x = lum(i);
% %     time = toc(tstart);
%     % Color the screen a random color
%     Screen('FillRect', window, [x x x]);
% %     DrawFormattedText(window,['Beep #1 \n' num2str(time)],'center','center',[0 0 1]);
%     % Flip to the screen
%     vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
% end
% PsychPortAudio('Stop',pahandle);
% 
% for i = 1:beepPauseLengthFrames
%     time = toc(tstart);
%     DrawFormattedText(window,['SILENCE #2'],'center','center',[0 0 1]);
%     Screen('Flip',window);
% end
% PsychPortAudio('Start',pahandle,repetitions,startCue,waitForDeviceStart);
% 
% vbl = Screen('Flip', window);
% for i = 1:beepLengthFrames
% %     DrawFormattedText(window,'BEEP #2','center','center',[0 0 1]);
% %     Screen('Flip',window);
%     x = lum(i);
% %     time = toc(tstart);
%     % Color the screen a random color
%     Screen('FillRect', window, [x x x]);
% %     DrawFormattedText(window,['Beep #2 \n' num2str(time)],'center','center',[0 0 1]);
% 
%     vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
% end
% PsychPortAudio('Stop',pahandle);
% 
% for i = 1:beepPauseLengthFrames
%     time = toc(tstart);
%     DrawFormattedText(window,['SILENCE #3'],'center','center',[0 0 1]);
%     Screen('Flip',window);
% end
% tstop = toc(tstart);
% 
% PsychPortAudio('Close',pahandle);
% sca;


%% Simple beep demo
% info = audiodevinfo;
% 
% InitializePsychSound(1);
% 
% nrchannels = 2;
% freq = 48000;
% 
% repetitions = 1;
% beepLengthSecs = 1;
% beepPauseTime = 1;
% 
% envelopeFreq = 1;
% 
% startCue = 0;
% 
% waitForDeviceStart = 1;
% 
% pahandle = PsychPortAudio('Open',-1,1,1,freq,nrchannels);
% 
% PsychPortAudio('Volume',pahandle,1);
% 
% % myBeep = MakeBeep(500,beepLengthSecs,freq);
% 
% % myBeep = zeros(beepLengthSecs*freq);
% myBeep = cos(440*2*pi/freq*(1:(beepLengthSecs*freq))).*(cos(envelopeFreq*2*pi/freq*(1:(beepLengthSecs*freq)))/2+0.5);
% 
% PsychPortAudio('FillBuffer',pahandle,[myBeep; myBeep]);
% 
% tstart = tic;
% 
% PsychPortAudio('Start',pahandle,repetitions,startCue,waitForDeviceStart);
% 
% [actualStartTime,endPositionSecs,xruns,estStopTime] = PsychPortAudio('Stop',pahandle,1,1);
% 
% startCue = estStopTime + beepPauseTime;
% 
% PsychPortAudio('Start',pahandle,repetitions,startCue,waitForDeviceStart);
% 
% PsychPortAudio('Stop',pahandle,1,1);
% 
% tstop = toc(tstart);
% tstart2 = tic;
% tstop2 = toc(tstart2);
% 
% PsychPortAudio('Close',pahandle);



%% Wait frames demo

% PsychDefaultSetup(2);
% 
% screens = Screen('Screens');
% 
% screenNumber = max(screens);
% 
% white = WhiteIndex(screenNumber);
% black = BlackIndex(screenNumber);
% grey = white/2;
% 
% % Screen('Preference', 'SkipSyncTests', 1);
% [window, windowRect] = PsychImaging('OpenWindow',screenNumber,grey);
% 
% ifi = Screen('GetFlipInterval', window);
% topPriorityLevel = MaxPriority(window);
% Priority(topPriorityLevel);
% 
% flipSecs = 1;
% waitframes = round(flipSecs / ifi);
% waitframes = 1;
% 
% vbl = Screen('Flip', window);
% 
% lum = cos(2*pi/60*(1:300))/2+0.5;
% 
% tstart = tic;
% for f = 1:300
%     x = lum(f);
%     % Color the screen a random color
%     Screen('FillRect', window, [x x x]);
% 
%     % Flip to the screen
%     vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
% 
% end
% tstop = toc(tstart);
% 
% sca;

%% Other stuff?
% rect = Screen('Rect', window);
% [screenXpixels, screenYpixels] = Screen('WindowSize', window);
% [xCenter, yCenter] = RectCenter(windowRect);
% ifi = Screen('GetFlipInterval', window);
% hertz = FrameRate(window);
% nominalHertz = Screen('NominalFrameRate', window);
% pixelSize = Screen('PixelSize', window);
% [width, height] = Screen('DisplaySize', screenNumber);
% maxLum = Screen('ColorRange', window);