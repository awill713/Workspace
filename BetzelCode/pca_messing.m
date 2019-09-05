% load fisheriris
load('D:\2Pdata\data\K070_20170803_2P_FRA_01.mat')

%%
% X = calcium.npilSubTraces';
X = dataM';
% X = double(imread('cameraman.tif'));
mu = mean(X);

[eigenvectors, scores,~,~,exp] = pca(X);

% Xrem = X-eigenvectors(:, 1)*(scores(:, 1)');

nComp = 2:size(scores,2);
Xhat = scores(:,nComp) * eigenvectors(:,nComp)';
Xhat = bsxfun(@plus, Xhat, mu);
figure
imagesc(X')
figure
imagesc(Xhat')
% X = X';
% Xhat = Xhat';


[bf,af] = butter(5,20/exptInfo.fr/2,'low'); % 10 is good for low pass
%         [sos,g] = tf2sos(bf,af);
stats.lowPassFilter.b = bf;
stats.lowPassFilter.a = af;
npst_f = zeros(size(calcium.npilSubTraces));
npst_z = npst_f;
for ii=1:size(calcium.npilSubTraces,1)
    %             npst_f(ii,:) = conv(double(calcium.npilSubTraces(ii,:)),sos,'same')*g;
    npst_f(ii,:) = filtfilt(bf,af,double(Xhat(:,ii)));
    npst_z(ii,:) = filtfilt(bf,af,double(X(:,ii)));
end

if iscell(stimInfo)
    stimInfo = stimInfo{1};
end


stimDur = stimInfo.tDur/1000;
ITI = stimInfo.ITI/1000;
preEv = floor(1*exptInfo.fr); % 1 second before each stim
postEv = ceil((stimDur+ITI)*exptInfo.fr); % x seconds post stim
if events.eventsOn(end)+postEv>length(npst_f)
    events.eventsOn = events.eventsOn(events.eventsOn+postEv<length(npst_f)==1);
    events.eventsOff = events.eventsOff(1:length(events.eventsOn)+1);
end
proc.preEv = preEv;
proc.postEv = postEv;
pcaRaster = makeCaRaster(npst_f,events.eventsOn,preEv,postEv,1);
caRaster =  makeCaRaster(npst_z,events.eventsOn,preEv,postEv,1);
% stats.rawRast = makeCaRaster(npst_f,events.eventsOn,preEv,postEv,0);

figure
plot(mean(caRaster(:,:,71)))
hold on
plot(mean(pcaRaster(:,:,71)))
figure
imagesc(caRaster(:,:,71))
figure
imagesc(pcaRaster(:,:,71))




