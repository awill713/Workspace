

file = '20181203AW056_dataStruct.mat';
folder = '/Volumes/AARON FILES/Two photon/2PM 004/AW056/AV WF 20181203/';

load([folder '/' file]);


mouse = exptInfo.mouse;
date = exptInfo.recDate;
frequencies = stimInfo.frequencies;
order = stimInfo.order;

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

offsetRange = -1:1:1;
entireMI = zeros(1,totalNeurons);
entireMItuned = zeros(1,length(offsetRange));
entireMIuntuned = zeros(1,length(offsetRange));
entireMIklds = zeros(1,length(offsetRange));
entireMIp = zeros(1,length(offsetRange));
entireMISet = zeros(totalNeurons,length(offsetRange));
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
    
    entireMISet(:,o) = entireMI;
    
    entireMItuned(o) = mean(entireMI(tunedNeurons));
    entireMIuntuned(o) = mean(entireMI(untunedNeurons));
    entireMIklds(o) = KLDivergence(entireMI(tunedNeurons),entireMI(untunedNeurons));
    [h, p] = ttest2(entireMI(tunedNeurons),entireMI(untunedNeurons));
    entireMIp(o) = p;
    offset
    
    shuffledContMItuned(o) = mean(shuffledContMI(tunedNeurons));
    shuffledContMIuntuned(o) = mean(shuffledContMI(untunedNeurons));
    
end