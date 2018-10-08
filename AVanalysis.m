

dimensions = [size(fluoresence.pixelFluor,1) size(fluoresence.pixelFluor,2)];
totalPixels = dimensions(1)*dimensions(2);

patchDim = [32 43];
range = -10:19;

patchNumber = totalPixels / (patchDim(1)*patchDim(2));

reshapedPixel = reshape(fluoresence.pixelFluor,[patchDim(1) totalPixels/patchDim(1) size(fluoresence.pixelFluor,3)]);

repeats = sum(events.order==2);
patchFluor = zeros(patchNumber,repeats,length(range));

lightOn = events.onset(events.order==2);

for p = 1:patchNumber
    p
    for r = 1:repeats
        timeOn = lightOn(r);
        trace = squeeze(mean(mean(squeeze(reshapedPixel(:,(patchDim(2)*(p-1)+1):patchDim(2)*p,range+timeOn)),2),1));
        baseline = mean(trace(1:10));
        patchFluor(p,r,:) = trace - baseline;
    end
end

figure;
imagesc(squeeze(patchFluor(:,1,:)));
xticks(1:5:30);xticklabels(range(1:5:end));
colorbar

    