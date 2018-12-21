
tunedNeurons = find([tuningStats.significant]);
untunedNeurons = find([tuningStats.significant] == 0);

figure;scatter(reshape(tuningCurves,[1720 1]),reshape(standardDeviations,[1720 1]));
xlabel('Tone response');
ylabel('Standard deviation');

figure;scatter(reshape(tuningCurves,[1720 1]),reshape(standardDeviations,[1720 1]).^2);
xlabel('Tone response');
ylabel('Variance');

figure;histogram([mutualInformation{:,2}],linspace(0,15,16));
legend('All neurons');
xlabel('Mutual information (bits)');
figure;histogram([mutualInformation{tunedNeurons,2}],linspace(0,15,16));
hold on;histogram([mutualInformation{untunedNeurons,2}],linspace(0,15,16));
legend('Tuned neurons','Untuned neurons');
xlabel('Mutual information (bits)');

[h p] = ttest2([mutualInformation{tunedNeurons,2}],[mutualInformation{untunedNeurons,2}]);
figure;hold on;
bar(1,mean([mutualInformation{tunedNeurons,2}]),'k');
bar(2,mean([mutualInformation{untunedNeurons,2}]),'r');

zeroedResponses = zeros(1,17200);
normalizedResponses = zeros(1,17200);
counter = 1;
for i = 1:totalNeurons
    for toneCount = 1:uniqueTones
        tempAvg = tuningCurves(i,toneCount);
        for r = 1:repeats
            zeroedResponses(1,counter) = mean(allResponses(i,toneCount,r,preToneFrames:(preToneFrames+framesForAverage))) - tempAvg;
            normalizedResponses(1,counter) = (mean(allResponses(i,toneCount,r,preToneFrames:(preToneFrames+framesForAverage))) - tempAvg)/(tempAvg*0.1532 + 0.3453);
            counter = counter + 1;
        end
    end
end
figure;histogram(zeroedResponses);
hold on;
gaus = normrnd(0,std(zeroedResponses),[1 17200]);
histogram(gaus);
legend('Actual zero-meaned responses','Generated gaussian');

figure;histogram(normalizedResponses);
hold on;
gaus = normrnd(0,std(normalizedResponses),[1 17200]);
histogram(gaus);
legend('Normalized responses','Generated gaussian');


% for n = 1:86
%     for i = 1:length(xpix{1,n})
%         tempImage(xpix{1,n}(i)+min(dat.ops.xrange),ypix{1,n}(i)+min(dat.ops.yrange)) = 1;
%     end
% end
% image = zeros(512,512,3);
% image(:,:,1) = tempImage;
% image(:,:,2) = tempImage;
% image(:,:,3) = tempImage;
cellNums = find([dat.stat.iscell]);
for i = 1:88
    n = cellNums(i);
    ipix{i} = dat.stat(n).ipix;
    ypix{i} = dat.stat(n).ypix;
    xpix{i} = dat.stat(n).xpix;
end

bFreq = zeros(1,88);
for n = 1:88
    bFreq(n) = find(tuningCurves(n,:)==max(tuningCurves(n,:)));
end

colorRange = jet;
coloredImage = zeros(512,512,3);
for i = 1:88
    if sum(ismember(tunedNeurons,i))
        freq = bFreq(i);
        temp = round(interp1(linspace(1,20,length(colorRange)),1:length(colorRange),freq));
        r = colorRange(temp,1);
        g = colorRange(temp,2);
        b = colorRange(temp,3);
    else
        r = 1;
        g = 1;
        b = 1;
    end
    
    xs = xpix{1,i};
    ys = ypix{1,i};
    
    for p = 1:length(xs)
        coloredImage(xs(p)+min(dat.ops.xrange),ys(p)+min(dat.ops.yrange),1) = r;
        coloredImage(xs(p)+min(dat.ops.xrange),ys(p)+min(dat.ops.yrange),2) = g;
        coloredImage(xs(p)+min(dat.ops.xrange),ys(p)+min(dat.ops.yrange),3) = b;
    end
end
figure;imshow(coloredImage);
colormap(jet);
cb = colorbar;
cb.Ticks = linspace(0,1,2);
cb.TickLabels = [4 70];
ylabel(cb,'Frequency (kHz)');

colorRange = jet;
coloredImage = ones(512,512,3).*0.5;
for i = 1:88
    if sum(ismember(tunedNeurons,i))
        freq = bFreq(i);
        temp = round(interp1(linspace(1,20,length(colorRange)),1:length(colorRange),freq));
        r = colorRange(temp,1);
        g = colorRange(temp,2);
        b = colorRange(temp,3);
    else
        r = 0;
        g = 0;
        b = 0;
    end
    
    
    
    xs = xpix{1,i};
    ys = ypix{1,i};
    
    for p = 1:length(xs)
        coloredImage(xs(p)+min(dat.ops.xrange),ys(p)+min(dat.ops.yrange),1) = r;
        coloredImage(xs(p)+min(dat.ops.xrange),ys(p)+min(dat.ops.yrange),2) = g;
        coloredImage(xs(p)+min(dat.ops.xrange),ys(p)+min(dat.ops.yrange),3) = b;
    end
end
figure;imshow(coloredImage);
colormap(jet);
cb = colorbar;
cb.Ticks = [0 5/20 10/20 15/20 1];
cb.TickLabels = [4 7.3 15 32 70];
ylabel(cb,'Frequency (kHz)');

