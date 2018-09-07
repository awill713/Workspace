
clear;
[file folder] = uigetfile('/Users/Aaron/Documents/');
if file==0
    return;
end
load([folder '/' file]);

mouse = exptInfo.mouse;
date = exptInfo.recDate;
frequencies = stimInfo.frequencies;

totalNeurons = length(calcium.n);
uniqueTones = length(stimInfo.index);

totalTones = length(stimInfo.order);
repeats = stimInfo.repeats;

frameRate = exptInfo.fr;
preToneTime = 1000; %ms
preToneFrames = round(frameRate*preToneTime/1000);
postToneTime = stimInfo.ITI; %ms
postToneFrames = ceil(frameRate*postToneTime/1000);

toneResponses = cell(totalNeurons,uniqueTones);
tuningCurves = zeros(totalNeurons,uniqueTones);
timeForAverage = 1000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);
     
for n = 1:totalNeurons
    for t = 1:totalTones
%         [n t]
        tone = stimInfo.order(t);
        firstFrame = events.eventsOn(t) - preToneFrames;
        lastFrame = events.eventsOn(t) + postToneFrames;
        
        pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(t)-1));
        pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(t)-1));
        trace = calcium.npilSubTraces(n,firstFrame:lastFrame)-pretoneAverage;
        trace = trace ./ pretoneSTD;

%         pretoneAverage = mean(calcium.rawTraces(n,firstFrame:events.eventsOn(t)-1));
%         pretoneSTD = std(calcium.rawTraces(n,firstFrame:events.eventsOn(t)-1));
%         trace = calcium.rawTraces(n,firstFrame:lastFrame)-pretoneAverage;
%         trace = trace ./ pretoneSTD;

        if isempty(toneResponses{n,tone}) %first instance of tone, so stored traces have yet to be initialized
            toneResponses{n,tone} = trace/repeats;
        else
            toneResponses{n,tone} = toneResponses{n,tone} + trace/repeats;
        end
    end
    
    for u = 1:uniqueTones
        tuningCurves(n,u) = mean(toneResponses{n,u}(preToneFrames:(preToneFrames+framesForAverage)));
    end

end

save([folder '/' mouse '_' date '_FRA'],'tuningCurves','toneResponses','frequencies');