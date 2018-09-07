function [ TE ] = TransferEntropy(xTrain,yTrain, offset )
%TRANSFERENTROPY Summary of this function goes here
%   Detailed explanation goes here

x = xTrain;
y = yTrain;
tau = offset;


tMinusTau = length(x)-tau;

jCounts = zeros(2);
probs = zeros(2);
for t = 1:tMinusTau
    t0 = y(t)+1;
    t1 = y(t+tau)+1;
    jCounts(t0,t1) = jCounts(t0,t1) + 1;
end
for c = 1:size(jCounts,2)
    for r = 1:size(jCounts,1)
        probs(r,c) = jCounts(r,c) / sum(jCounts(r,:));
    end
end
ijCounts = zeros(2,2,2);
ijProbs = zeros(2,2,2);
for t = 1:tMinusTau
    i0 = x(t)+1;
    j0 = y(t)+1;
    j1 = y(t+tau)+1;
    ijCounts(i0,j0,j1) = ijCounts(i0,j0,j1) + 1;
end
for z = 1:size(ijCounts,3)
    for c = 1:size(ijCounts,2)
        for r = 1:size(ijCounts,1)
            ijProbs(r,c,z) = ijCounts(r,c,z) / sum(ijCounts(r,c,:));
        end
    end
end

TE = 0;

for gI = 1:size(ijProbs,1)
    for gJ = 1:size(ijProbs,2)
        for r = 1:size(ijProbs,3)
            gI;
            gJ;
            r;
            overallProb = ijCounts(gI,gJ,r)/sum(sum(sum(ijCounts)));
            if overallProb == 0
                dTE = 0;
            else
                dTE = overallProb * log2(ijProbs(gI,gJ,r)/probs(gJ,r));
            end
            TE = TE + dTE;
        end
    end
end

end

