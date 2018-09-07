

stimTrain = randi([1 9],1,1000);
mid = 5;
rSTD = 1;
respTrain = zeros(1,length(stimTrain));

for i = 1:length(stimTrain)
    respTrain(i) = 1/(rSTD*sqrt(2*pi)) * exp(-0.5*((stimTrain(i)-mid)/rSTD)^2) + randn*0.05;
end

uniqueStimuli = length(unique(stimTrain));
stimuliValues = unique(stimTrain);



gaussRgivenS = zeros(uniqueStimuli,3); %stimulus state, response mean, response std

for i = 1:length(stimuliValues)
    tempStimValue = stimuliValues(i);
    tempPoints = zeros(1,sum(stimTrain==i));
    
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

figure;
errorbar(gaussRgivenS(:,1),gaussRgivenS(:,2),gaussRgivenS(:,3));
