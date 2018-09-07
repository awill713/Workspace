function [ TE averageShuffled significance sigRatio sigDifference ] = TransferEntropyCalcium( x, y, bins, delay )
%TRANSFERENTROPYCALCIUM Summary of this function goes here
%   Detailed explanation goes here

repeats = 100;

xTrain = x;
yTrain = y;
binNumber = bins;
offset = delay;

xBinEdges = linspace(min(xTrain),max(xTrain),binNumber+1);
binnedX = discretize(xTrain,xBinEdges);

yBinEdges = linspace(min(yTrain),max(yTrain),binNumber+1);
binnedY = discretize(yTrain,yBinEdges);

TE = calculateTE(binnedX, binnedY, binNumber, offset);


shuffledTE = zeros(1,repeats);
for repeat = 1:repeats
    shuffledX = binnedX(randperm(length(binnedX)));
    deltaAverage = calculateTE(shuffledX,binnedY,binNumber,offset);
    shuffledTE(repeat) = deltaAverage;
    repeat
end
averageShuffled = mean(shuffledTE);
% std(shuffledTE)
threshold = averageShuffled + std(shuffledTE);
if TE > threshold
    significance = 1;
else
    significance = 0;
end
sigRatio = TE / threshold;
sigDifference = TE - threshold;

end

function [transferEntropy] = calculateTE(train1,train2,binCount,tau)
tMinusTau = length(train2) - tau;

yCounts = zeros(binCount);
for t = 1:(length(train2)-1)
    t0 = train2(t);
    t1 = train2(t+1);
    yCounts(t0,t1) = yCounts(t0,t1)+1;
end
yProbs = zeros(binCount);
for r = 1:size(yCounts,1)
    rowSum = sum(yCounts(r,:));
    for c = 1:size(yCounts,2)
        yProbs(r,c) = yCounts(r,c) / rowSum;
    end
end

xyCounts = zeros(binCount,binCount,binCount);
for t = 1:tMinusTau
    x0 = train1(t);
    y0 = train2(t+tau-1);
    y1 = train2(t+tau);
    xyCounts(x0,y0,y1) = xyCounts(x0,y0,y1)+1;
end
xyProbs = zeros(binCount,binCount,binCount);
for i = 1:size(xyCounts,1)
    for j = 1:size(xyCounts,2)
        surfaceSum = sum(xyCounts(i,j,:));
        for k = 1:size(xyCounts,3)
            xyProbs(i,j,k) = xyCounts(i,j,k)/surfaceSum;
        end
    end
end

transferEntropy = 0;

for teX = 1:size(xyProbs,1)
    for teY = 1:size(xyProbs,2)
        for teR = 1:size(xyProbs,3)
            jointProb = xyCounts(teX,teY,teR) / sum(sum(sum(xyCounts)));
            conditionalProb = xyProbs(teX,teY,teR);
            responseProb = yProbs(teY,teR);
            if jointProb == 0
                dTE = 0;
            else
                dTE = jointProb * log2(conditionalProb / responseProb);
            end
            transferEntropy = transferEntropy + dTE;
        end
    end
end
end

