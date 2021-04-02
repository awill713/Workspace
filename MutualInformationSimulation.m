
clear
lambda = 0:.01:1;
trials = 10;
orientations = 12;
repeats = 100;

info = zeros(2,length(lambda)); %mean and std

for lam = 1:length(lambda)
    lam
    tempInfo = zeros(1,repeats);
    for r = 1:repeats
    spikes = poissrnd(lambda(lam),[orientations trials]);
    
    peak = max(spikes(:));
    valley = min(spikes(:));
    subHist = zeros(orientations,10);
    negEntropy = zeros(1,orientations);
    for o = 1:orientations
        subHist(o,:) = histcounts(spikes(o,:),linspace(valley,peak+1,11),'Normalization','probability');
        subEntropy = subHist(o,:).*log2(subHist(o,:)); subEntropy(isnan(subEntropy))=0;
        
        negEntropy(o) = sum(subEntropy);
    end
    totalHist = sum(subHist)./orientations;
    
    entropy = totalHist.*log2(totalHist); entropy(isnan(entropy))=0;
    totalEntropy = -sum(entropy);
    
    
    tempInfo(r) = totalEntropy - (-sum(1/orientations * negEntropy));
    end
    info(1,lam) = mean(tempInfo);
    info(2,lam) = std(tempInfo);
end

figure;hold on;
m = info(1,:);
s = info(2,:);
hA = area(lambda,[m-s; 2*s]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(lambda,m,'Color',[0 0 0]);
