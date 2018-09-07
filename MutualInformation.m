function [ s2rInfo s2rInfoTemp ] = MutualInformation( s, r, c, d )
%MUTUALINFORMATION Summary of this function goes here
%   Detailed explanation goes here

stimTrain = s;
respTrain = r;
binCount = c;

uniqueStimuli = length(unique(stimTrain));
stimuliValues = unique(stimTrain);

responseBinEdges = linspace(min(respTrain),max(respTrain),binCount+1);
binnedResponse = discretize(respTrain,responseBinEdges);

responseProbability = zeros(1,binCount);
for i = 1:binCount
    responseProbability(i) = sum(binnedResponse==i)/length(binnedResponse);
end


%%%Figure 3

countRgivenS = zeros(uniqueStimuli,binCount);
probRgivenS = zeros(uniqueStimuli,binCount);

gaussRgivenS = zeros(uniqueStimuli,4); %stimulus state, response mean, response std, normalization scalar

for i = 1:length(stimuliValues)
    tempStimValue = stimuliValues(i);
    tempPoints = zeros(1,sum(stimTrain==tempStimValue));
    
    stimIndex = find(stimTrain==tempStimValue);
    for j = 1:length(stimIndex)
        tempPoints(j) = respTrain(stimIndex(j));
    end
    
    tempMean = mean(tempPoints);
    tempSTD = std(tempPoints);
    gaussRgivenS(i,1) = tempStimValue;
    gaussRgivenS(i,2) = tempMean;
    gaussRgivenS(i,3) = tempSTD;
end

% figure;
% errorbar(gaussRgivenS(:,1),gaussRgivenS(:,2),gaussRgivenS(:,3));


%calculate overall probability of response
lowest = min(respTrain) - std(respTrain)/2;
highest = max(respTrain) + std(respTrain)/2;
increment = std(respTrain)/1000;
totalRespProb = zeros(2,length(lowest:increment:highest)); %response value and probability of that response value

totalRespProb(1,:) = lowest:increment:highest;
for r = 1:length(totalRespProb)
    tempR = 0;
    
    for s = 1:uniqueStimuli
        avg = gaussRgivenS(s,2);
        sDev = gaussRgivenS(s,3);
        dR = 1/(sDev*sqrt(2*pi)) * exp(-0.5*((totalRespProb(1,r)-avg)/sDev)^2);
        
        tempR = tempR + dR;
    end
    totalRespProb(2,r) = tempR / uniqueStimuli;
end

%scale the probabilities for each stimulus such that they add to 1
% figure;
for s = 1:uniqueStimuli
    avg = gaussRgivenS(s,2);
    sDev = gaussRgivenS(s,3);
    rx = lowest:increment:highest;
    yx = zeros(1,length(rx));
    for rr = 1:length(rx)
        yx(rr) = 1/(sDev*sqrt(2*pi)) * exp(-0.5*((rx(rr)-avg)/sDev)^2);
    end
    probScalar = 1/sum(yx);
    gaussRgivenS(s,4) = probScalar;
    yx = yx*probScalar;
%     plot(rx,yx);
%     hold on;
    sum(yx);
end
% legend('1','2','3','4','5','6','7','8','9','10');

totalRespProb(2,:) = totalRespProb(2,:)/sum(totalRespProb(2,:));
% figure;
% plot(totalRespProb(1,:),totalRespProb(2,:));

% totalRespProb(1,:)
% totalRespProb(2,:)

%calculate MI for each stimulus
s2rInfoTemp = zeros(1,uniqueStimuli);

for s = 1:length(stimuliValues)
    tempS2Rinfo = 0;
    for r = 1:length(totalRespProb)
        theMean = gaussRgivenS(s,2);
        theSTD = gaussRgivenS(s,3);
        scalar = gaussRgivenS(s,4);
        theRgivenS = scalar * 1/(theSTD*sqrt(2*pi)) * exp(-0.5*((totalRespProb(1,r)-theMean)/theSTD)^2);
        
        dMI = theRgivenS * log2(theRgivenS / totalRespProb(2,r));
        tempS2Rinfo = tempS2Rinfo + dMI;
    end
    s2rInfoTemp(s) = tempS2Rinfo;
end

% figure;
% plot(s2rInfoTemp);

for i = 1:length(stimTrain)
    stim = stimTrain(i);
    resp = binnedResponse(i);
    
    countRgivenS(stim,resp) = countRgivenS(stim,resp) + 1;
end

for s = 1:size(countRgivenS,1)
    for r = 1:size(countRgivenS,2)
        probRgivenS(s,r) = countRgivenS(s,r) / sum(countRgivenS(s,:));
    end
end

s2rInfo = zeros(1,uniqueStimuli);

for s = 1:length(s2rInfo)
    tempInfo = 0;
    for r = 1:length(responseProbability)
        if probRgivenS(s,r)==0
            dInfo = 0;
        else
            dInfo = probRgivenS(s,r) * log2(probRgivenS(s,r)/responseProbability(r));
        end
        tempInfo = tempInfo + dInfo;
    end
    s2rInfo(s) = tempInfo;
end

MI = sum(s2rInfo);

end

