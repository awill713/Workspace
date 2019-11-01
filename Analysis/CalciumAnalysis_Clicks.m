
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
eachTrial = zeros(totalNeurons,repeats,totalFrames);
timeForAverage = 1000; %ms;
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
    
    eachTrial(n,:,:) = tempResp;
    
end

responsiveNeurons = find([tuningStats.significant]==1);
length(responsiveNeurons)
length(find([tuningStats.responseSign]==1))
length(find([tuningStats.responseSign]==-1))

up = responsiveNeurons(find([tuningStats.responseSign]==1));
down = responsiveNeurons(find([tuningStats.responseSign]==-1));

f1 = figure;hold on;
plot(mean(responses(up,:)));
plot(mean(responses(down,:)));
legend({['n = ' num2str(length(up))], ['n = ' num2str(length(down))]});
title('Sound responsive neurons');
xticks(1:round(frameRate):totalFrames);
xticklabels({'-1','0','1','2','3','4','5','6','7'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');

save([folder '/' mouse '_' date '_clickResponses'],'responses','eachTrial','tuningStats','responsiveNeurons');
saveas(f1,[folder '/' mouse '_' date '_clickResponsesFigure.fig']);