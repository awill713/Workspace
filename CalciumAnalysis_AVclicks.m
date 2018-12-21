
clear;
[file folder] = uigetfile('/Users/Aaron/Documents/'); %looking for *processedData.mat
if file==0
    return;
end
load([folder '/' file]);

mouse = exptInfo.mouse;
date = exptInfo.recDate;

totalNeurons = length(calcium.n);

% totalTones = length(stimInfo.order);
repeats = stimInfo.repeats;

frameRate = exptInfo.fr;
preToneTime = 1000; %ms
preToneFrames = round(frameRate*preToneTime/1000);
postToneTime = stimInfo.interClickInterval*1000; %ms
postToneFrames = ceil(frameRate*postToneTime/1000);
totalFrames = postToneFrames + preToneFrames + 1;

responses = zeros(totalNeurons,3,totalFrames);
timeForAverage = 1000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

stimTypes = unique(stimInfo.order);

for n = 1:totalNeurons
    
%     respForStats = zeros(repeats,uniqueTones);
%     respForStats2 = zeros(repeats,uniqueTones);
    tempResp = zeros(3,repeats,totalFrames);
    
    for type = 1:length(stimTypes)
        eventsOfType = find(stimInfo.order == type);
        for r = 1:repeats
            ev = eventsOfType(r);
            firstFrame = events.eventsOn(ev) - preToneFrames;
            lastFrame = events.eventsOn(ev) + postToneFrames;
            
            pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(ev)-1));
            pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(ev)-1));
            
            trace = calcium.npilSubTraces(n,firstFrame:lastFrame) - pretoneAverage;
            trace = trace ./ pretoneSTD;
            
            tempResp(type,r,:) = trace;
        end
    end
    
%     tempResp = tempResp ./ repeats;
    
    responses(n,:,:) = squeeze(mean(tempResp,2));
    
    soundPre = mean(tempResp(1,:,1:preToneFrames),3);
    lightPre = mean(tempResp(2,:,1:preToneFrames),3);
    bothPre = mean(tempResp(3,:,1:preToneFrames),3);
    
    soundPost = mean(tempResp(1,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    lightPost = mean(tempResp(2,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    bothPost = mean(tempResp(3,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    
    pre = mean(tempResp(:,1:preToneFrames),2);
    post = mean(tempResp(:,preToneFrames+1:preToneFrames+1+framesForAverage),2);
    
    
    [tuningStats(n).soundSig tuningStats(n).soundpValue] = ttest(soundPost);
    [tuningStats(n).lightSig tuningStats(n).lightpValue] = ttest(lightPost);
    [tuningStats(n).bothSig tuningStats(n).bothpValue] = ttest(bothPost);
    
    [tuningStats(n).significant tuningStats(n).pValue] = ttest(post);
    if tuningStats(n).soundSig
        tuningStats(n).soundResponseSign = sign(mean(soundPost)-mean(soundPre));
    end
    if tuningStats(n).lightSig
        tuningStats(n).lightResponseSign = sign(mean(lightPost)-mean(lightPre));
    end
    if tuningStats(n).bothSig
        tuningStats(n).bothResponseSign = sign(mean(bothPost)-mean(bothPre));
    end
    
    tuningStats(n).AVinteractpValue = anova1([soundPost' lightPost' bothPost'],[],'off');
    tuningStats(n).AVinteractSig = tuningStats(n).AVinteractpValue < 0.05;
    
    
end

soundResponsiveNeurons = find([tuningStats.soundSig]==1);
lightResponsiveNeurons = find([tuningStats.lightSig]==1);
bothResponsiveNeurons = find([tuningStats.bothSig]==1);
% length(responsiveNeurons)
% length(find([tuningStats.responseSign]==1))
% length(find([tuningStats.responseSign]==-1))


save([folder '/' mouse '_' date '_AVclickResponses'],'responses','tuningStats','soundResponsiveNeurons','lightResponsiveNeurons','bothResponsiveNeurons');