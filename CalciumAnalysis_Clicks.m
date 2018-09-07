
clear;
[file folder] = uigetfile('/Users/Aaron/Documents/');
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

responses = zeros(totalNeurons,totalFrames);
timeForAverage = 2000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

for n = 1:totalNeurons
    
%     respForStats = zeros(repeats,uniqueTones);
%     respForStats2 = zeros(repeats,uniqueTones);
    tempResp = zeros(repeats,totalFrames);
    
    for r = 1:repeats
        
        firstFrame = events.eventsOn(r) - preToneFrames;
        lastFrame = events.eventsOn(r) + postToneFrames;
        
        pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(r)-1));
        pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(r)-1));
        
        trace = calcium.npilSubTraces(n,firstFrame:lastFrame) - pretoneAverage;
        trace = trace ./ pretoneSTD;
        
        tempResp(r,:) = trace;
    end
    
%     tempResp = tempResp ./ repeats;
    
    responses(n,:) = mean(tempResp);
    
    pre = mean(tempResp(:,1:preToneFrames),2);
    post = mean(tempResp(:,preToneFrames+1:preToneFrames+1+framesForAverage),2);
    
    [tuningStats(n).significant tuningStats(n).pValue] = ttest(post);
    if tuningStats(n).significant
        tuningStats(n).responseSign = sign(mean(post)-mean(pre));
    end
    
end

responsiveNeurons = find([tuningStats.significant]==1);
length(responsiveNeurons)
length(find([tuningStats.responseSign]==1))
length(find([tuningStats.responseSign]==-1))

% meanTrace = zeros(1,totalFrames);
% posTrace = zeros(1,totalFrames);
% posCount = 0;
% negTrace = zeros(1,totalFrames);
% negCount = 0;
% for nn = 1:length(responsiveNeurons)
%     neuron = responsiveNeurons(nn);
%     if tuningStats(neuron).responseSign > 0
%         posTrace = posTrace + responses(neuron,:);
%         posCount = posCount+1;
%     else
%         negTrace = negTrace + responses(neuron,:);
%         negCount = negCount+1;
%     end
%     meanTrace = meanTrace + responses(neuron,:);
% end
% figure;
% subplot(1,3,1);
% posTrace = posTrace./posCount;
% plot(posTrace);
% subplot(1,3,2);
% negTrace = negTrace./negCount;
% plot(negTrace);
% subplot(1,3,3);
% meanTrace = meanTrace./length(responsiveNeurons);
% plot(meanTrace);

save([folder '/' mouse '_' date '_clickResponses'],'responses','tuningStats','responsiveNeurons');