
clear;

experiment = 'EP001';
mouseID = 'AW093';
session = 'Session1';
date = '20191029';
stimPath = 'C:\Users\Aaron\Documents\MATLAB\Workspace\Analysis\EphyStimInfo\20190623_multiAmpTones_opto_ONE_REPEAT_70dB_400k_005_stimInfo';
analysisWindow = [-50 200];
quantWindow = [0 30]; %ms after stimulus onset
rasterHeatMap = 'jet';

% dataFolder = fullfile('D:\KiloSort\',mouseID,session,date,'SpikeMat');
dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*multAmpTones_opto*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,date,'OptoMultiAmpToneResponses');
newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'OptoMultiAmpToneResponses');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(stimPath);
frequencies = stimInfo.frequencies;
amplitudes = stimInfo.intensities;
uniqueFreq = length(frequencies);
uniqueAmp = length(amplitudes);
uniqueEvents = uniqueFreq*uniqueAmp*2;

totalUnits = length(dataFiles);

frScalar = 1000/(quantWindow(2)-quantWindow(1));
map = colormap(rasterHeatMap);close;

%ERROR IN STIMULUS GENERATION LED TO LASER ON AT INDEX 100 (and laser=2 at
%index 200)
stimInfo.index(100,4)=0;
stimInfo.index(200,4)=1;

for n = 84:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    tuningRaster = cell(uniqueFreq,uniqueAmp,2,3); %spike times (1), trial number (y-axis) (2), color (3)
    tuningMat = zeros(uniqueFreq,uniqueAmp,2);
    for f = 1:uniqueFreq
        freq = frequencies(f);
        freqIndices = find(stimInfo.index(:,2)==freq);
        color = map(round(f/length(frequencies)*length(map)),:);
        
        for a = 1:uniqueAmp
            amp = amplitudes(a);
            ampIndices = find(stimInfo.index(:,3)==amp);
            
            for laser = 1:2
                laserCondition = laser-1;
                laserIndices = find(stimInfo.index(:,4)==laserCondition);
                
                event1Index = intersect(intersect(freqIndices,ampIndices),laserIndices);
                event1 = find(stimInfo.order==event1Index);
                
                eventsOfInterest = event1:uniqueEvents:length(eventTimes);
                
                spikeCount = zeros(length(eventsOfInterest),1);
                rasterX = [];
                rasterY = [];
                tempY = (a-1)*length(eventTimes)/uniqueAmp + (f-1)*length(eventTimes)/(uniqueEvents/2) +1;
                for e = 1:length(eventsOfInterest)
                    time = eventTimes(eventsOfInterest(e));
                    alignedSpikes = (spikeTimes-time)/1000; %temporal resolution of cheetah is in microseconds, converting to milliseconds
                    trialSpikes = alignedSpikes(find(alignedSpikes>analysisWindow(1) & alignedSpikes<analysisWindow(end)));
                    
                    quant = length(find(trialSpikes>quantWindow(1) & trialSpikes<quantWindow(2)));
                    spikeCount(e,1) = quant;
                    
                    tempY = tempY+1;
                    rasterX = [rasterX trialSpikes];
                    rasterY = [rasterY repmat(tempY,[1 length(trialSpikes)])];
                end
                
                tuningRaster{f,a,laser,1} = rasterX;
                tuningRaster{f,a,laser,2} = rasterY;
                tuningRaster{f,a,laser,3} = color;
                tuningMat(f,a,laser) = mean(spikeCount)*frScalar;
            end
        end
    end
    
    noLaserVector = reshape(tuningMat(:,:,1),[uniqueFreq*uniqueAmp 1]);
    yesLaserVector = reshape(tuningMat(:,:,2),[uniqueFreq*uniqueAmp 1]);
    fit = polyfit(noLaserVector,yesLaserVector,1);
    
    unitData(n).raster = tuningRaster;
    unitData(n).freqAmpOptoTuning = tuningMat;
    unitData(n).optoEffect = fit;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    
    low = min(min(min(tuningMat)));
    high = max(max(max(tuningMat)));
    
    f1 = figure;
    set(f1,'Position',[100 50 800 600])
    subplot(2,3,1); colormap 'jet';
    for r = 1:size(tuningRaster,1)
        for col = 1:size(tuningRaster,2)
            scatter(tuningRaster{r,col,1,1},tuningRaster{r,col,1,2},[],tuningRaster{r,col,1,3},'.');
            hold on;
        end
    end
    for a = 1:uniqueAmp-1
        line([analysisWindow(1) analysisWindow(end)],[a*uniqueEvents/2/uniqueAmp a*uniqueEvents/2/uniqueAmp],'LineStyle',':');
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial number and type');
    c = colorbar;
    set(c,'YTick',[0 1]);
    set(c,'YTickLabel',[frequencies(1,1) frequencies(1,end)])
    ylabel(c,'Frequency');
    title('No laser');
    %     title(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6))]);
    
    subplot(2,3,4);
    imagesc(flip(tuningMat(:,:,1)'));
    xlist = [1 5:5:20];
    xticks(xlist);
    xticklabels(round(frequencies(xlist)/1000,1))
    xlabel('Frequency (kHz)');
    ylist = 1:2:uniqueAmp;
    yticks(ylist);
    yticklabels(flip(amplitudes(ylist)));
    ylabel('SPL (dB)');
    cb = colorbar;
    caxis([low high]);
    ylabel(cb,'Firing rate (Hz)');
    title('No laser');
    
    subplot(2,3,2); colormap 'jet';
    for r = 1:size(tuningRaster,1)
        for col = 1:size(tuningRaster,2)
            scatter(tuningRaster{r,col,2,1},tuningRaster{r,col,2,2},[],tuningRaster{r,col,2,3},'.');
            hold on;
        end
    end
    for a = 1:uniqueAmp-1
        line([analysisWindow(1) analysisWindow(end)],[a*uniqueEvents/2/uniqueAmp a*uniqueEvents/2/uniqueAmp],'LineStyle',':');
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial number and type');
    c = colorbar;
    set(c,'YTick',[0 1]);
    set(c,'YTickLabel',[frequencies(1,1) frequencies(1,end)])
    ylabel(c,'Frequency');
    title('With laser');
    
    subplot(2,3,5);
    imagesc(flip(tuningMat(:,:,2)'));
    xlist = [1 5:5:20];
    xticks(xlist);
    xticklabels(round(frequencies(xlist)/1000,1))
    xlabel('Frequency (kHz)');
    ylist = 1:2:uniqueAmp;
    yticks(ylist);
    yticklabels(flip(amplitudes(ylist)));
    ylabel('SPL (dB)');
    cb = colorbar;
    caxis([low high]);
    ylabel(cb,'Firing rate (Hz)');
    title('With laser');
    
%     subplot(2,3,3);
%     imagesc(flip((tuningMat(:,:,2)-tuningMat(:,:,1))'));
%     xlist = [1 5:5:20];
%     xticks(xlist);
%     xticklabels(round(frequencies(xlist)/1000,1))
%     xlabel('Frequency (kHz)');
%     ylist = 1:2:uniqueAmp;
%     yticks(ylist);
%     yticklabels(flip(amplitudes(ylist)));
%     ylabel('SPL (dB)');
%     cb = colorbar;
%     ylabel(cb,'Firing rate (Hz)');
%     title('Difference');
    
    subplot(2,3,3);
    for spl = 1:size(tuningMat,2)
        semilogx(frequencies/1000,tuningMat(:,spl,2) - tuningMat(:,spl,1));
        hold on;
    end
    line([1 10^ceil(log10(frequencies(end)/1000))],[0 0],'LineStyle','--');
    xlabel('Frequency (kHz)');
    ylabel('\Delta FR (Hz)');
    legend(strsplit(num2str(amplitudes)));
    
%     set(f1,'Position',[150 150 1200 600]);
%     
    subplot(2,3,6);
    scatter(noLaserVector,yesLaserVector);
    unity = refline(1,0); unity.Color = 'k'; unity.LineStyle = '--';
    fitLine = refline(fit(1),fit(2)); fitLine.Color = 'b';
    xlabel('Tone evoked FR, without laser (Hz)');
    ylabel('Tone evoked FR, with laser (Hz)');
    title(['Slope = ' num2str(fit(1)) ', y-int = ' num2str(fit(2))]);
    
    suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6))]);
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.jpg']));
    close(f1);
end

f2 = figure;
subplot(1,2,1);
slopes = [unitData.optoEffect]; slopes = slopes(1:2:end);
intercepts = [unitData.optoEffect]; intercepts = intercepts(2:2:end);
histogram(slopes,20);
title('Slopes');
subplot(1,2,2);
histogram(intercepts,20);
title('y-intercepts');
saveas(f2,fullfile(newDir,'Slope and intercept summary'));

analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.folder = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow; %ms
analysisParams.quantificationWindow = quantWindow;
save(fullfile(newDir,'OptoMultiAmpToneData.mat'),'unitData','analysisParams');
