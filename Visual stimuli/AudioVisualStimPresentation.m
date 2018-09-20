
sca;
close all;
clear all;

%% Create auditory stimulus and events

sampleFreq = 48000;

freqMin = 220;
freqMax = 880;
uniqueFreq = 10;
repeats = 10;

toneDuration = 1; %seconds
interToneInterval = 5; %seconds between tone end and tone onset

freqRange = 10.^linspace(log10(freqMin),log10(freqMax),uniqueFreq);

order = [];
for r = 1:repeats
    order = [order randperm(uniqueFreq)];
end

%% Beep text demo

%Sound setup
InitializePsychSound(1);

nrchannels = 2;
freq = 48000;

soundFreq = 440;
envelopeFreq = 1;

repetitions = 1;
beepLengthSecs = 5;
beepPauseTime = 2;
startCue = 0;
preStimSilence = 10;

waitForDeviceStart = 1;
pahandle = PsychPortAudio('Open',[],3,1,freq,[2 1]);
PsychPortAudio('Volume',pahandle,1);
% myBeep = MakeBeep(500,beepLengthSecs,freq);
myBeep = cos(soundFreq*2*pi/freq*(1:(beepLengthSecs*freq))).*(cos(envelopeFreq*2*pi/freq*(1:(beepLengthSecs*freq)))/2+0.5);
silence = zeros(1,beepPauseTime*freq);
preSilence = zeros(1,preStimSilence*freq);
stimulus = [preSilence myBeep silence myBeep silence];

eventsOut = zeros(1,length(stimulus));
eventsOut(1) = 5;

PsychPortAudio('FillBuffer',pahandle,[stimulus;stimulus]);

a = GetSecs;
x= PsychPortAudio('GetStatus',pahandle)
PsychPortAudio('GetAudioData',pahandle,100);
PsychPortAudio('Start',pahandle,repetitions,startCue,waitForDeviceStart);

PsychPortAudio('Stop',pahandle,1,1);

framesOut = PsychPortAudio('GetAudioData',pahandle);

b = GetSecs;
PsychPortAudio('Close',pahandle);
sca;
b-a

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