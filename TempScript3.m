

%% Load and process data
% clear;
file = 'AW017_20180519_2P_FRA_processedData.mat';
folder = '/Volumes/AARON FILES/Two photon/2PM 001/AW017/2P 20180519/Pure tones/';
% [file folder] = uigetfile('/Users/Aaron/Documents/');
% if file==0
%     return;
% end
load([folder '/' file]);

mouse = exptInfo.mouse;
date = exptInfo.recDate;
frequencies = stimInfo.frequencies;
order = stimInfo.order;

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

toneResponses = cell(totalNeurons,uniqueTones);
standardDeviations = zeros(totalNeurons,uniqueTones);
tuningCurves = zeros(totalNeurons,uniqueTones);
mutualInformation = cell(totalNeurons,2);
timeForAverage = 2000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

shuffledToneResponses = cell(totalNeurons,uniqueTones);
shuffledTuningCurves = zeros(totalNeurons,uniqueTones);
shuffledSTD = zeros(totalNeurons,uniqueTones);
shuffledMI = cell(totalNeurons,2);

timeIterations = totalFrames - preToneFrames;
variedDurResp = zeros(totalNeurons,uniqueTones,timeIterations);
variedDurSTD = zeros(totalNeurons,uniqueTones,timeIterations);
variedDurMI = zeros(totalNeurons,timeIterations);


allDurResp = zeros(totalNeurons,uniqueTones,151,151);
allDurSTD = zeros(totalNeurons,uniqueTones,151,151);
allDurMI = zeros(totalNeurons,151,151);

pastTime = 2000; %ms
pastTimeFrames = ceil(frameRate*pastTime/1000);
totalExpFrames = length(calcium.npilSubTraces);
processedFluor = zeros(totalNeurons,totalExpFrames);
toneFluors = zeros(totalNeurons,uniqueTones,repeats,totalFrames);
mutualInformationFluor = cell(totalNeurons,2);
variedDurRespFluor = zeros(totalNeurons,uniqueTones,timeIterations);
variedDurSTDFluor = zeros(totalNeurons,uniqueTones,timeIterations);
variedDurMIFluor = zeros(totalNeurons,timeIterations);
avgToneFluor = cell(totalNeurons,uniqueTones);
tuningCurvesFluor = zeros(totalNeurons,uniqueTones);
standardDeviationsFluor = zeros(totalNeurons,uniqueTones);


for n = 1:totalNeurons
    
    shuffledTrace(n,:) = calcium.npilSubTraces(n,randperm(length(calcium.npilSubTraces)));
%     shuffledTrace = rand([1, length(calcium.npilSubTraces)]);
    
    respForStats = zeros(repeats,uniqueTones);
    respForStatsFluor = zeros(repeats,uniqueTones);
    
    for u = 1:uniqueTones
        
        uTones = find(order==u);
        
        for t = 1:length(uTones)
            uT = uTones(t);
            firstFrame = events.eventsOn(uT) - preToneFrames;
            lastFrame = events.eventsOn(uT) + postToneFrames;
            
            pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(uT)-1));
            pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(uT)-1));
        
            trace = calcium.npilSubTraces(n,firstFrame:lastFrame)-pretoneAverage;
            trace = trace ./ pretoneSTD;

            allResponses(n,u,t,:) = trace;
            
            pretoneAverage = mean(shuffledTrace(n,firstFrame:events.eventsOn(uT)-1));
            pretoneSTD = std(shuffledTrace(n,firstFrame:events.eventsOn(uT)-1));
            
            trace = shuffledTrace(n,firstFrame:lastFrame)-pretoneAverage;
            trace = trace ./ pretoneSTD;
            
            shuffledResponses(n,u,t,:) = trace;
        end

        toneResponses{n,u} = squeeze(mean(allResponses(n,u,:,:)))';
        tuningCurves(n,u) = mean(toneResponses{n,u}(preToneFrames:(preToneFrames+framesForAverage)));
        standardDeviations(n,u) = std(squeeze(mean(allResponses(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4)));
        standardErrors(n,u) = standardDeviations(n,u)./sqrt(repeats);
        respForStats(:,u) = squeeze(mean(allResponses(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4));
        
        shuffledToneResponses{n,u} = squeeze(mean(shuffledResponses(n,u,:,:)))';
        shuffledTuningCurves(n,u) = mean(shuffledToneResponses{n,u}(preToneFrames:(preToneFrames+framesForAverage)));
        shuffledSTD(n,u) = std(squeeze(mean(shuffledResponses(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4)));
        
        for td = 1:timeIterations
            variedDurResp(n,u,td) = mean(toneResponses{n,u}(preToneFrames:(preToneFrames+td)));
            variedDurSTD(n,u,td) = std(squeeze(mean(allResponses(n,u,:,preToneFrames:(preToneFrames+td)),4)));
        end
        
%         for startAvg = 1:timeIterations
%             timeLeft = timeIterations - startAvg;
%             for endAvg = 1:timeLeft
%                 allDurResp(n,u,startAvg,startAvg+endAvg) = mean(toneResponses{n,u}((preToneFrames+startAvg):(preToneFrames+startAvg+endAvg)));
%                 allDurSTD(n,u,startAvg,startAvg+endAvg) = std(squeeze(mean(allResponses(n,u,:,(preToneFrames+startAvg):(preToneFrames+startAvg+endAvg)),4)));
%             end
%         end
    end
 
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
    
    [curve, value] = MutualInformation2(tuningCurves(n,:),standardDeviations(n,:));
    mutualInformation{n,1} = curve;
    mutualInformation{n,2} = value;
    
    [curve, value] = MutualInformation2(shuffledTuningCurves(n,:),shuffledSTD(n,:));
    shuffledMI{n,1} = curve;
    shuffledMI{n,2} = value;
    
    for td = 1:timeIterations
        [c, v] = MutualInformation2(variedDurResp(n,:,td),variedDurSTD(n,:,td));
        variedDurMI(n,td) = v;
    end
    
%     for startAvg = 1:timeIterations
%         timeLeft = timeIterations - startAvg;
%         for endAvg = 1:timeLeft
%             [c v] = MutualInformation2(allDurResp(n,:,startAvg,startAvg+endAvg),allDurSTD(n,:,startAvg,startAvg+endAvg));
%             allDurMI(n,startAvg,startAvg+endAvg) = v;
%         end
%     end
    
    for f = pastTimeFrames+1 : totalExpFrames
        pastF = mean(calcium.npilSubTraces(n,(f-pastTimeFrames):(f-1)));
        pastSTD = std(calcium.npilSubTraces(n,(f-pastTimeFrames):(f-1)));
        processedFluor(n,f) = (calcium.npilSubTraces(n,f) - pastF) / pastSTD;
    end
    for u = 1:uniqueTones
        
        uTones = find(order==u);
        
        for t = 1:length(uTones)
            uT = uTones(t);
            firstFrame = events.eventsOn(uT) - preToneFrames;
            lastFrame = events.eventsOn(uT) + postToneFrames;
            
            trace = processedFluor(n,firstFrame:lastFrame);

            toneFluors(n,u,t,:) = trace;
            raw(n,u,t,:) = calcium.npilSubTraces(n,firstFrame:lastFrame);
        end

        avgToneFluor{n,u} = squeeze(mean(toneFluors(n,u,:,:)))';
        
        tuningCurvesFluor(n,u) = mean(avgToneFluor{n,u}(preToneFrames:(preToneFrames+framesForAverage)));
        standardDeviationsFluor(n,u) = std(squeeze(mean(toneFluors(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4)));
        standardErrorsFluor(n,u) = standardDeviationsFluor(n,u)./sqrt(repeats);
        respForStatsFluor(:,u) = squeeze(mean(toneFluors(n,u,:,preToneFrames:(preToneFrames+framesForAverage)),4));
        
        for td = 1:timeIterations
            variedDurRespFluor(n,u,td) = mean(avgToneFluor{n,u}(preToneFrames:(preToneFrames+td)));
            variedDurSTDFluor(n,u,td) = std(squeeze(mean(toneFluors(n,u,:,preToneFrames:(preToneFrames+td)),4)));
        end
        
%         for startAvg = 1:timeIterations
%             timeLeft = timeIterations - startAvg;
%             for endAvg = 1:timeLeft
%                 allDurResp(n,u,startAvg,startAvg+endAvg) = mean(toneResponses{n,u}((preToneFrames+startAvg):(preToneFrames+startAvg+endAvg)));
%                 allDurSTD(n,u,startAvg,startAvg+endAvg) = std(squeeze(mean(allResponses(n,u,:,(preToneFrames+startAvg):(preToneFrames+startAvg+endAvg)),4)));
%             end
%         end
    end
    tuningStatsFluor(n).pValue = anova1(respForStatsFluor,[],'off');
    tuningStatsFluor(n).significant = tuningStatsFluor(n).pValue<0.05;
    
    if tuningStatsFluor(n).significant
        sig = zeros(1,repeats);
        sig1Samp = zeros(1,repeats);
        
        for uu = 1:uniqueTones
            pre = squeeze(mean(toneFluors(n,uu,:,1:preToneFrames-1),4));
            post = squeeze(mean(toneFluors(n,uu,:,preToneFrames:(preToneFrames+framesForAverage)),4));
            
            sig(uu) = ttest(pre,post);
            sig1Samp(uu) = ttest(post);
        end
        
        tuningStatsFluor(n).sigTones = find(sig~=0);
        tuningStatsFluor(n).sigTones1Samp = find(sig1Samp~=0);
    end
    
    [curve, value] = MutualInformation2(tuningCurvesFluor(n,:),standardDeviationsFluor(n,:));
    mutualInformationFluor{n,1} = curve;
    mutualInformationFluor{n,2} = value;
    
    for td = 1:timeIterations
        [c, v] = MutualInformation2(variedDurRespFluor(n,:,td),variedDurSTDFluor(n,:,td));
        variedDurMIFluor(n,td) = v;
    end
    n
end

tunedNeurons = find([tuningStats.significant]);
untunedNeurons = find([tuningStats.significant] == 0);

%% Display basic MI data
figure;
subplot(2,1,1); hold on;
histogram([mutualInformation{tunedNeurons,2}],linspace(0,10,21));
histogram([mutualInformation{untunedNeurons,2}],linspace(0,10,21));
histogram([shuffledMI{tunedNeurons,2}],linspace(0,10,21));
histogram([shuffledMI{untunedNeurons,2}],linspace(0,10,21));
legend('Tuned neurons','Untuned neurons','Shuffled tuned','Shuffled untuned');
xlabel('Mutual information (bits)');

subplot(2,1,2); hold on;
[h, p] = ttest2([mutualInformation{tunedNeurons,2}],[mutualInformation{untunedNeurons,2}]);
m1 = mean([mutualInformation{tunedNeurons,2}]); s1 = std([mutualInformation{tunedNeurons,2}]);
m2 = mean([mutualInformation{untunedNeurons,2}]); s2 = std([mutualInformation{untunedNeurons,2}]);
m3 = mean([shuffledMI{tunedNeurons,2}]); s3 = std([shuffledMI{tunedNeurons,2}]);
m4 = mean([shuffledMI{untunedNeurons,2}]); s4 = std([shuffledMI{untunedNeurons,2}]);
bar(1,mean([mutualInformation{tunedNeurons,2}]),'b');
bar(2,mean([mutualInformation{untunedNeurons,2}]),'r');
bar(3,mean([shuffledMI{tunedNeurons,2}]),'k');
bar(4,mean([shuffledMI{untunedNeurons,2}]),'g');
errorbar([1 2 3 4],[m1 m2 m3 m4], [s1 s2 s3 s4],'.');
title(['p-value = ' num2str(p)]);


%% MI based on varied framesForAverage
tunedMIs = mean(variedDurMI(tunedNeurons,:),1);
untunedMIs = mean(variedDurMI(untunedNeurons,:),1);
figure; 
subplot(2,2,1); hold on;
plot(tunedMIs);
plot(untunedMIs);
xlabel('Frames for average');
ylabel('Mutual information (bits)');
legend('Tuned neurons','Untuned neurons');
title('Mutual information');

subplot(2,2,2);
plot(tunedMIs - untunedMIs);
xlabel('Frames for average');
ylabel('\Delta mutual information (bits)');
title('MI difference');
subplot(2,2,3);
plot(tunedMIs./untunedMIs);
xlabel('Frames for average');
ylabel('MI_{tuned} / MI_{untuned}');
title('MI ratio');

klds = zeros(1,timeIterations);
pvalues = zeros(1,timeIterations);
for td = 1:timeIterations
    klds(td) = KLDivergence(variedDurMI(tunedNeurons,td),variedDurMI(untunedNeurons,td));
    [h, p] = ttest2(variedDurMI(tunedNeurons,td),variedDurMI(untunedNeurons,td));
    pvalues(td) = p;
end
subplot(2,2,4);
plot(klds);
yyaxis right;
plot(log10(pvalues));
legend('KL Divergence','p-value')
xlabel('Frames for average');
title('KL Divergence');

%% MI based on varied start and end points of averaging
% klds3D = zeros(151,151);
% for startAvg = 1:timeIterations
%     timeLeft = timeIterations - startAvg;
%     for endAvg = 1:timeLeft
%         tempTuned = allDurResp(tunedNeurons,:,startAvg,startAvg+endAvg);
%         tempUntuned = allDurResp(unTunedNeurons,:,startAvg,startAvg+endAvg);
%         tempKLD = KLDivergence(tempTuned,tempUntuned);
%         klds3D(startAvg,startAvg+endAvg) = tempKLD;
%     end
% end
% figure;
% surface(klds3D);

%% Basic fluor analysis
figure; 
subplot(2,2,1);
plot(toneResponses{40,7})
title('Normal analysis, n40t7');
subplot(2,2,2);
plot(avgToneFluor{40,7})
title('Sliding baseline, n40t7');

subplot(2,2,3);
plot(toneResponses{73,9})
title('Normal analysis, n73t9');
subplot(2,2,4);
plot(avgToneFluor{73,9})
title('Sliding baseline, n73t9');

%% Mutual information from entire trace, equidistant bins
bins = uniqueTones+1;
stimulusMaster = zeros(1,totalExpFrames);
for e = 1:length(events.eventsOn)
    stimulusMaster(events.eventsOn(e):events.eventsOff(e+1)) = order(e);
end
% stimulus = randi(21,[1 32000])-1;
% range = linspace(-0.5,max(stimulus)+0.5,bins+1);
% stimHist = histcounts(stimulus,range);
% stimProb = stimHist ./ sum(stimHist);
% stimEntropy = sum(-stimProb .* log2(stimProb));

offsetRange = -25:1:25;
entireMI = zeros(1,totalNeurons);
entireMItuned = zeros(1,length(offsetRange));
entireMIuntuned = zeros(1,length(offsetRange));
entireMIklds = zeros(1,length(offsetRange));
entireMIp = zeros(1,length(offsetRange));
peakMIandOffset = zeros(2,totalNeurons);

shuffledContMI = zeros(1,totalNeurons);
shuffledContMItuned = zeros(1,length(offsetRange));
shuffledContMIuntuned = zeros(1,length(offsetRange));


for o = 1:length(offsetRange)
    
    offset = offsetRange(o);
    
    stimulus = stimulusMaster;
    if offset>0
        stimulus((end-offset+1):end) = [];
    elseif offset<0
        stimulus(1:-offset) = [];
    end
    range = linspace(-0.5,max(stimulus)+0.5,bins+1);
    stimHist = histcounts(stimulus,range);
    stimProb = stimHist ./ sum(stimHist);
    stimEntropy = sum(-stimProb .* log2(stimProb));
    
    entireMI = zeros(1,totalNeurons);
    
    for n = 1:totalNeurons
        responseTrace = calcium.npilSubTraces(n,:);
        shuff = shuffledTrace(n,:);
%         responseTrace = processedFluor(n,:);
%         responseTrace = rand([1 32000]);
%         responseTrace = shuffledTrace(n,:);
        
        minFluor = min(responseTrace);
        maxFluor = max(responseTrace);
        binnedFluor = discretize(responseTrace,linspace(minFluor,maxFluor,bins+1));
        binnedShuff = discretize(shuff,linspace(minFluor,maxFluor,bins+1));
        
        entropyGivenFluor = zeros(1,length(unique(binnedFluor)));
        fluorProb = zeros(1,length(unique(binnedFluor)));
        
        entropyGivenFluorShuff = zeros(1,length(unique(binnedFluor)));
        shuffProb = zeros(1,length(unique(binnedFluor)));
        
        for fluor = min(binnedFluor):max(binnedFluor)

            indices = find(binnedFluor == fluor);
            if offset>0
                indices(indices<offset+1) = [];
            elseif offset<0
                indices(indices>(length(responseTrace)+offset)) = []; 
            end
            stimGivenFluor = stimulusMaster(indices-offset);
            tempHist = histcounts(stimGivenFluor,range);
            tempProb = tempHist ./ sum(tempHist);
            tempProb(isnan(tempProb)) = 0;
            entropyGivenFluor(fluor) = sum(-tempProb .* log2(tempProb+eps));
            
            fluorProb(fluor) = sum(binnedFluor == fluor) / length(binnedFluor);
            
            indices = find(binnedShuff == fluor);
            if offset>0
                indices(indices<offset+1) = [];
            elseif offset<0
                indices(indices>(length(responseTrace)+offset)) = []; 
            end
            stimGivenFluor = stimulusMaster(indices-offset);
            tempHist = histcounts(stimGivenFluor,range);
            tempProb = tempHist ./ sum(tempHist);
            tempProb(isnan(tempProb)) = 0;
            entropyGivenFluorShuff(fluor) = sum(-tempProb .* log2(tempProb+eps));
            
            shuffProb(fluor) = sum(binnedShuff == fluor) / length(binnedShuff);
        end
        
        tempMI = stimEntropy - sum(fluorProb.* entropyGivenFluor);
        tempMIShuff = stimEntropy - sum(shuffProb.* entropyGivenFluorShuff);
        
        entireMI(n) = tempMI;
        shuffledContMI(n) = tempMIShuff;
        if entireMI(n) > peakMIandOffset(1,n)
            peakMIandOffset(1,n) = entireMI(n);
            peakMIandOffset(2,n) = offset;
        end
    end
    
    entireMItuned(o) = mean(entireMI(tunedNeurons));
    entireMIuntuned(o) = mean(entireMI(untunedNeurons));
    entireMIklds(o) = KLDivergence(entireMI(tunedNeurons),entireMI(untunedNeurons));
    [h, p] = ttest2(entireMI(tunedNeurons),entireMI(untunedNeurons));
    entireMIp(o) = p;
    offset
    
    shuffledContMItuned(o) = mean(shuffledContMI(tunedNeurons));
    shuffledContMIuntuned(o) = mean(shuffledContMI(untunedNeurons));
    
end
figure;
subplot(2,1,1);hold on;
plot(offsetRange,entireMItuned);
plot(offsetRange,entireMIuntuned);
plot(offsetRange,shuffledContMItuned);
plot(offsetRange,shuffledContMIuntuned);
subplot(2,1,2); hold on;
plot(offsetRange,entireMItuned./entireMIuntuned);
plot(offsetRange,ones(1,length(offsetRange)));
figure;
plot(offsetRange,entireMIklds);
yyaxis right;
plot(offsetRange,log(entireMIp));
legend('KL Divergence','p-value')
xlabel('Offset');
title('Continuous MI');
figure; hold on;
histogram(peakMIandOffset(2,tunedNeurons),-300:15:300);
histogram(peakMIandOffset(2,untunedNeurons),-300:15:300);
% figure; hold on;
% [h, p] = ttest2(entireMI(tunedNeurons),entireMI(untunedNeurons));
% bar(1,mean(entireMI(tunedNeurons)),'b');
% bar(2,mean(entireMI(untunedNeurons)),'r');
% title(['p-value = ' num2str(p)]);

%% Mutual information for entire trace, varying temporal summation
bins = uniqueTones+1;
stimulusMaster = zeros(1,totalExpFrames);
for e = 1:length(events.eventsOn)
    stimulusMaster(events.eventsOn(e):events.eventsOff(e+1)) = order(e);
end
% stimulus = randi(21,[1 32000])-1;
% range = linspace(-0.5,max(stimulus)+0.5,bins+1);
% stimHist = histcounts(stimulus,range);
% stimProb = stimHist ./ sum(stimHist);
% stimEntropy = sum(-stimProb .* log2(stimProb));

framesBack = 0:25:1500;
temporalMI = zeros(1,totalNeurons);
temporalMItuned = zeros(1,length(framesBack));
temporalMIuntuned = zeros(1,length(framesBack));
temporalMIklds = zeros(1,length(framesBack));
temporalMIp = zeros(1,length(framesBack));
temporalPeakMIandOffset = zeros(2,totalNeurons);

for fb = 1:length(framesBack)
    
    retroFrames = framesBack(fb);
    
    stimulus = stimulusMaster(retroFrames+1:end);
    range = linspace(-0.5,max(stimulus)+0.5,bins+1);
    stimHist = histcounts(stimulus,range);
    stimProb = stimHist ./ sum(stimHist);
    stimEntropy = sum(-stimProb .* log2(stimProb));
    
    temporalMI = zeros(1,totalNeurons);
    
    for n = 1:totalNeurons
        responseTrace = zeros(1,totalExpFrames-retroFrames);
        index = 1;
        for rf = retroFrames+1:totalExpFrames
            responseTrace(index) = mean(calcium.npilSubTraces(n,rf-retroFrames:rf));
            index = index+1;
        end
        
        minFluor = min(responseTrace);
        maxFluor = max(responseTrace);
        binnedFluor = discretize(responseTrace,linspace(minFluor,maxFluor,bins+1));
        
        entropyGivenFluor = zeros(1,length(unique(binnedFluor)));
        fluorProb = zeros(1,length(unique(binnedFluor)));
        
        for fluor = min(binnedFluor):max(binnedFluor)

            indices = find(binnedFluor == fluor);
            stimGivenFluor = stimulusMaster(indices);
            tempHist = histcounts(stimGivenFluor,range);
            tempProb = tempHist ./ sum(tempHist);
            tempProb(isnan(tempProb)) = 0;
            entropyGivenFluor(fluor) = sum(-tempProb .* log2(tempProb+eps));
            
            fluorProb(fluor) = sum(binnedFluor == fluor) / length(binnedFluor);
        end
        
        tempMI = stimEntropy - sum(fluorProb.* entropyGivenFluor);
        temporalMI(n) = tempMI;
        if temporalMI(n) > temporalPeakMIandOffset(1,n)
            temporalPeakMIandOffset(1,n) = temporalMI(n);
            temporalPeakMIandOffset(2,n) = offset;
        end
        n
    end
    
    temporalMItuned(fb) = mean(temporalMI(tunedNeurons));
    temporalMIuntuned(fb) = mean(temporalMI(untunedNeurons));
    temporalMIklds(fb) = KLDivergence(temporalMI(tunedNeurons),temporalMI(untunedNeurons));
    [h, p] = ttest2(temporalMI(tunedNeurons),temporalMI(untunedNeurons));
    temporalMIp(fb) = p;
    retroFrames

end

figure; hold on;
plot(framesBack, temporalMItuned);
plot(framesBack, temporalMIuntuned);