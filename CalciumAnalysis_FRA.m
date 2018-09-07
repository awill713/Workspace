
clear;
[file folder] = uigetfile('/Users/Aaron/Documents/');
if file==0
    return;
end
load([folder '/' file]);

mouse = exptInfo.mouse;
date = exptInfo.recDate;
frequencies = stimInfo.frequencies;
% attenuations = stimInfo.attenuations;
order = stimInfo.order;
% index = stimInfo.index;

totalNeurons = length(calcium.n);
uniqueTones = stimInfo.uniqueTones;

totalTones = length(stimInfo.order);
repeats = stimInfo.repeats;

frameRate = exptInfo.fr;
preToneTime = 1000; %ms
preToneFrames = round(frameRate*preToneTime/1000);
postToneTime = stimInfo.ITI; %ms
postToneFrames = ceil(frameRate*postToneTime/1000);
totalFrames = postToneFrames + preToneFrames + 1;

% toneResponses = cell(totalNeurons,length(attenuations));
toneResponses = cell(totalNeurons,uniqueTones);
% tuningCurves = zeros(totalNeurons,length(attenuations),uniqueTones);
standardDeviations = zeros(totalNeurons,uniqueTones);
tuningCurves = zeros(totalNeurons,uniqueTones);
mutualInformation = cell(totalNeurons,2);
timeForAverage = 2000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

% for n = 1:totalNeurons
%     
%     for a = 1:length(attenuations)
%         
%         toneResponses{n,a} = zeros(uniqueTones,totalFrames);
%         
%         theResponses = zeros(uniqueTones,repeats,totalFrames);
%         respForStats = zeros(repeats,uniqueTones);
%         
%         aTones = find(index(:,3)==attenuations(a));
%         
%         for u = 1:length(aTones)
%             
%             uTones = find(order==aTones(u));
% 
%             for r = 1:length(uTones)
%                 uT = uTones(r);
%                 firstFrame = events.eventsOn(uT) - preToneFrames;
%                 lastFrame = events.eventsOn(uT) + postToneFrames;
%                 
%                 pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(uT)-1));
%                 pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(uT)-1));
%                 
%                 trace = calcium.npilSubTraces(n,firstFrame:lastFrame)-pretoneAverage;
%                 trace = trace ./ pretoneSTD;
% 
%                 theResponses(u,r,:) = trace;
%             end
%             
%             toneResponses{n,a}(u,:) = squeeze(mean(theResponses(u,:,:)))';
%             
%             tuningCurves(n,a,u) = mean(toneResponses{n,a}(u,preToneFrames:(preToneFrames+framesForAverage)));
% 
%             respForStats(:,u) = squeeze(mean(theResponses(u,:,preToneFrames:(preToneFrames+framesForAverage)),3));
%         end
%         
%         tuningStats(n).pValue(a) = anova1(respForStats,[],'off');
%         tuningStats(n).significant(a) = tuningStats(n).pValue(a)<0.05;
%         
%         if tuningStats(n).significant(a)
%             sig = zeros(1,repeats);
%             sig1Samp = zeros(1,repeats);
%             
%             for uu = 1:uniqueTones
%                 pre = squeeze(mean(theResponses(uu,:,1:preToneFrames-1),3));
%                 post = squeeze(mean(theResponses(uu,:,preToneFrames:(preToneFrames+framesForAverage)),3));
%                 
%                 
%                 sig(uu) = ttest(post);
%             end
%             
%             tuningStats(n).sigTones{a} = find(sig~=0);
%         end
%     end
% end



for n = 1:totalNeurons
    
    respForStats = zeros(repeats,uniqueTones);
    respForStats2 = zeros(repeats,uniqueTones);
    
    for u = 1:uniqueTones
        
        uTones = find(order==u);
        
%         tempResponses = zeros(repeats,totalFrames);
        
        for t = 1:length(uTones)
            uT = uTones(t);
            firstFrame = events.eventsOn(uT) - preToneFrames;
            lastFrame = events.eventsOn(uT) + postToneFrames;
            
            pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(uT)-1));
            pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(uT)-1));
        
            trace = calcium.npilSubTraces(n,firstFrame:lastFrame)-pretoneAverage;
            trace = trace ./ pretoneSTD;
            
%             tempResponses(t,:) = trace;
            
            allResponses(n,u,t,:) = trace;
        end
        
%         toneResponses{n,u} = mean(tempResponses);
%         length(toneResponses{n,u})
        toneResponses{n,u} = squeeze(mean(allResponses(n,u,:,:)))';
%         isequal(mean(tempResponses),toneResponses{n,u})
        
        tuningCurves(n,u) = mean(toneResponses{n,u}(preToneFrames:(preToneFrames+framesForAverage)));
%         standardDeviations(n,u) = std(toneResponses{n,u}(preToneFrames:(preToneFrames+framesForAverage)));
        standardDeviations(n,u) = std(squeeze(mean(allResponses(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4)));
        standardErrors(n,u) = standardDeviations(n,u)./sqrt(repeats);
%         respForStats(:,u) = mean(tempResponses(:,preToneFrames:(preToneFrames+framesForAverage)),2);
        respForStats(:,u) = squeeze(mean(allResponses(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4));
    end
%     
%     tuningStats(n).pValue = anova1(respForStats,[],'off'); %1 way anova without groups and display off
%     tuningStats(n).significant = tuningStats(n).pValue<0.05;
%     
    tuningStats(n).pValue = anova1(respForStats,[],'off');
    tuningStats(n).significant = tuningStats(n).pValue<0.05;
    
    if tuningStats(n).significant
        sig = zeros(1,repeats);
        sig1Samp = zeros(1,repeats);
        
        for uu = 1:uniqueTones
            pre = squeeze(mean(allResponses(n,uu,:,1:preToneFrames-1),4));
            post = squeeze(mean(allResponses(n,uu,:,preToneFrames:(preToneFrames+framesForAverage)),4));
            
            
            sig(uu) = ttest(pre,post);
            sig1Samp(uu) = ttest(post);
        end
        
        tuningStats(n).sigTones = find(sig~=0);
        tuningStats(n).sigTones1Samp = find(sig1Samp~=0);
    end
    
    [curve value] = MutualInformation2(tuningCurves(n,:),standardDeviations(n,:));
    mutualInformation{n,1} = curve;
    mutualInformation{n,2} = value;
end

% anovaMatrix = zeros(totalNeurons*repeats,uniqueTones);
% for n = 1:totalNeurons
%     for uu = 1:uniqueTones
%         uList = find(order==uu);
%         for r = 1:length(uList)
%             temp = mean(allResponses(n,uu,r,preToneFrames:(preToneFrames+framesForAverage)));
%             anovaMatrix((n-1)*repeats+r,uu) = temp;
%         end
%     end
% end


% tempMat = reshape([tuningStats.significant],[length(attenuations) totalNeurons])';
% for a = 1:length(attenuations)
%     tunedNeurons{a} = find(tempMat(:,a)==1);
% end

tunedNeurons = find([tuningStats.significant]);


% save([folder '/' mouse '_' date '_FRA'],'tuningCurves','toneResponses','frequencies','attenuations','tuningStats','tunedNeurons');