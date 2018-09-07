function [ xcorrOutput maxOffset averageShuffled sig ratioSig diffSig ] = CrossCorrCalcium( x, y, offset )
%CROSSCORRCALCIUM Summary of this function goes here
%   Detailed explanation goes here

xTrain = x;
yTrain = y;

xcorrArray = xcorr(xTrain,yTrain,offset,'coeff');
xcorrOutput = max(xcorrArray);
maxOffset = find(xcorrArray == xcorrOutput);

repeats = 100;
shuffledXCorr = zeros(1,repeats);
for repeat = 1:repeats
    shuffledX = xTrain(randperm(length(xTrain)));
    shuffledY = yTrain(randperm(length(yTrain)));
    shuffledXCorr(repeat) = max(xcorr(shuffledX,shuffledY,offset,'coeff'));
end
averageShuffled = mean(shuffledXCorr);
threshold = averageShuffled + std(shuffledXCorr);
if xcorrOutput > threshold
    sig = 1;
else
    sig = 0;
end
ratioSig = xcorrOutput / averageShuffled;
diffSig = xcorrOutput - averageShuffled;

end

