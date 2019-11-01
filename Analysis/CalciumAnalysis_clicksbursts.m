
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
postToneTime = stimInfo.ISI; %ms
postToneFrames = ceil(frameRate*postToneTime/1000);
totalFrames = postToneFrames + preToneFrames + 1;

responses = zeros(totalNeurons,2,totalFrames);
timeForAverage = 1000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

stimTypes = unique(stimInfo.order);

for n = 1:totalNeurons
    
%     respForStats = zeros(repeats,uniqueTones);
%     respForStats2 = zeros(repeats,uniqueTones);
    tempResp = zeros(2,repeats,totalFrames);
    
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
    
    clickPre = mean(tempResp(1,:,1:preToneFrames),3);
    burstPre = mean(tempResp(2,:,1:preToneFrames),3);
    
    clickPost = mean(tempResp(1,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    burstPost = mean(tempResp(2,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    
%     pre = mean(tempResp(:,1:preToneFrames),2);
%     post = mean(tempResp(:,preToneFrames+1:preToneFrames+1+framesForAverage),2);
    
    
    [tuningStats(n).clickSig tuningStats(n).clickpValue] = ttest(clickPost);
    [tuningStats(n).burstSig tuningStats(n).burstpValue] = ttest(burstPost);
    
%     [tuningStats(n).significant tuningStats(n).pValue] = ttest(post);
    if tuningStats(n).clickSig
        tuningStats(n).clickResponseSign = sign(mean(clickPost)-mean(clickPre));
    end
    if tuningStats(n).burstSig
        tuningStats(n).burstResponseSign = sign(mean(burstPost)-mean(burstPre));
    end
    
end

clickResponsiveNeurons = find([tuningStats.clickSig]==1);
burstResponsiveNeurons = find([tuningStats.burstSig]==1);
% length(responsiveNeurons)
% length(find([tuningStats.responseSign]==1))
% length(find([tuningStats.responseSign]==-1))


save([folder '/' mouse '_' date '_clickburstResponses'],'responses','tuningStats','clickResponsiveNeurons','burstResponsiveNeurons');

barResponses = zeros(totalNeurons,3);

for n = 1:totalNeurons
    for e = 1:2
        tri = zeros(1,20);
         
        barResponses(n,e) = squeeze(mean(responses(n,e,31:60),3));
    end
end
            
    


