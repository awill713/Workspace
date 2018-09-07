
sca;
close all;
clear all;

PsychDefaultSetup(2);

screens = Screen('Screens');

screenNumber = max(screens);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white/2;

Screen('Preference', 'SkipSyncTests', 1);
[window, windowRect] = PsychImaging('OpenWindow',screenNumber,white);

ifi = Screen('GetFlipInterval', window);
topPriorityLevel = MaxPriority(window);
numSecs = 1;
numFrames = round(numSecs / ifi);
flipSecs = 1;
waitframes = round(flipSecs / ifi);

vbl = Screen('Flip', window);

% Run until a key is pressed
while ~KbCheck

    % Color the screen a random color
    Screen('FillRect', window, rand(1, 3));

    % Flip to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

end

Priority(0);

sca;

% rect = Screen('Rect', window);
% [screenXpixels, screenYpixels] = Screen('WindowSize', window);
% [xCenter, yCenter] = RectCenter(windowRect);
% ifi = Screen('GetFlipInterval', window);
% hertz = FrameRate(window);
% nominalHertz = Screen('NominalFrameRate', window);
% pixelSize = Screen('PixelSize', window);
% [width, height] = Screen('DisplaySize', screenNumber);
% maxLum = Screen('ColorRange', window);