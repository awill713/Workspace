

channelCount = 5;
samples = 1000;
referenceChannel = 1;

%% Create raw signals of random noise and plot
raw = zeros(channelCount,samples);
for c = 1:channelCount
    raw(c,:) = rand(1,samples);
end
figure;hold on;
plot(raw');
title('Raw channels')

%% reference raw channels to "referenceChannel" and plot
referenced = zeros(channelCount,samples);
for c = 1:channelCount
    referenced(c,:) = raw(c,:) - raw(referenceChannel,:);
end
figure;hold on;
plot(referenced');
title(['Referenced to channel ' num2str(referenceChannel)]);

%% REreference the RAW based on its commmon average and plot
comAvgRefRaw = zeros(channelCount,samples);
comMedRefRaw = zeros(channelCount,samples);
commonAverageRaw = mean(raw,1);
commonMedianRaw = median(raw,1);
for c = 1:channelCount
    comAvgRefRaw(c,:) = raw(c,:) - commonAverageRaw;
    comMedRefRaw(c,:) = raw(c,:) - commonMedianRaw;
end
figure;hold on;
plot(comAvgRefRaw');
title('Raw, subtract the common average')

%% REreference the REFERENCED based on its common average and plot
comAvgRefRef = zeros(channelCount,samples);
comMedRefRef = zeros(channelCount,samples);
commonAverageRef = mean(referenced,1);
commonMedianRef = median(referenced,1);
for c = 1:channelCount
    comAvgRefRef(c,:) = referenced(c,:) - commonAverageRef;
    comMedRefRef(c,:) = referenced(c,:) - commonMedianRef;
end
figure;hold on;
plot(comAvgRefRef');
title(['Referenced to channel ' num2str(referenceChannel) ' then subtract the common average'])

%% look at correlation between common average (or common median) reference with and without first subtracting the referenceChannel

noSubtractRefChanMean = reshape(comAvgRefRaw,[1 samples*channelCount]);
yesSubtractRefChanMean = reshape(comAvgRefRef,[1 samples*channelCount]);
figure;scatter(noSubtractRefChanMean,yesSubtractRefChanMean)
xlabel(['Did NOT reference to channel ' num2str(referenceChannel) ' first']);
ylabel(['Did reference to channel ' num2str(referenceChannel) ' first']);
title('Common average subtracted');

noSubtractRefChanMedian = reshape(comMedRefRaw,[1 samples*channelCount]);
yesSubtractRefChanMedian = reshape(comMedRefRef,[1 samples*channelCount]);
figure;scatter(noSubtractRefChanMedian,yesSubtractRefChanMedian)
xlabel(['Did NOT reference to channel ' num2str(referenceChannel) ' first']);
ylabel(['Did reference to channel ' num2str(referenceChannel) ' first']);
title('Common median sutracted');

figure;scatter(noSubtractRefChanMean,noSubtractRefChanMedian)
xlabel('Rereferenced to common average');
ylabel('Rereferenced to common median');