
clear;

mouseID = 'AW000';
session = 'Session1';
folder = '2019-06-19_18-56-06';
stimPath = 'D:\Code\TuningCurve\TuningCurveSpkrLeft1-70k';
binWindows = -50:1:200; %ms
analysisWindow = [10 30]; %ms after stimulus onset
rasterHeatMap = 'jet';

dataFolder = fullfile('D:\KiloSort\',mouseID,session,folder,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'TuningCurveSpkr*'));

newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'TuningCurves','Figures');
mkdir(newDir);

load(stimPath);
frequencies = STIM.freqOrder;
frequencies(2,:) = 1:length(frequencies); %order presented
frequencies = sortrows(frequencies')';
uniqueTones = length(frequencies);

totalUnits = length(dataFiles);

binEdges = binWindows*1000; %temporal resolution of cheetah acquisition is microseconds
analysisEdges = analysisWindow*1000;
frScalar = 1000/(analysisWindow(2)-analysisWindow(1));
map = colormap(rasterHeatMap);

for n = 1:totalUnits
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    tuningRaster = cell(uniqueTones,3); %spike times, trial number (y-axis), color
    tuningMat = zeros(uniqueTones,length(binEdges)-1);
    tuningCurve = zeros(uniqueTones,1);
    tempY = 0;
    for f = 1:length(frequencies)
        freqIndex = frequencies(2,f);
        freqEvents = freqIndex:uniqueTones:length(eventTimes);
        
        trials = zeros(length(freqEvents),length(binEdges)-1);
        quantities = zeros(length(freqEvents),1);
        rasterX = [];
        rasterY = [];
        color = map(round(f/length(frequencies)*length(map)),:);
        for e = 1:length(freqEvents)
            time = eventTimes(freqEvents(e));
            counts = histcounts(spikeTimes,binEdges+time);
            quant = length(find(spikeTimes>(time+analysisEdges(1)) & spikeTimes<(time+analysisEdges(2))));
            trials(e,:) = counts;
            quantities(e,1) = quant;
            
            tempY = tempY+1;
            alignedSpikes = spikeTimes-time;
            trialSpikes = alignedSpikes(find(alignedSpikes>binEdges(1) & alignedSpikes<binEdges(end)))/1000;
            rasterX = [rasterX trialSpikes];
            rasterY = [rasterY repmat(tempY,[1 length(trialSpikes)])];
        end
        tuningRaster{f,1} = rasterX;
        tuningRaster{f,2} = rasterY;
        tuningRaster{f,3} = color;
        tuningMat(f,:) = mean(trials)*frScalar;
        tuningCurve(f) = mean(quantities)*frScalar;
    end
    
    unitData(n).raster = tuningRaster;
    unitData(n).PSTH = tuningMat;
    unitData(n).tuningCurve = tuningCurve;
    unitData(n).type = nData.CellInfo(6);
    
    f1 = figure;
    subplot(3,1,1); colormap 'jet';
    for r = 1:size(tuningRaster,1)
        scatter(tuningRaster{r,1},tuningRaster{r,2},[],tuningRaster{r,3},'.');
        hold on;
    end
    xlim([binWindows(1) binWindows(end)]);
    xlabel('Time (ms)');
    ylabel('Trial number and type');
    c = colorbar;
    set(c,'YTick',[0 1]);
    set(c,'YTickLabel',[frequencies(1,1) frequencies(1,end)])
    ylabel(c,'Frequency');
    title(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6))]);
    
    subplot(3,1,2);
    imagesc(flip(tuningMat,1));
    xlist = 1:50:250;
    xticks(xlist);
    xticklabels(binWindows(xlist))
    xlabel('Time (ms)');
    ylist = [1 10:10:50];
    yticks(ylist);
    yticklabels(flip(frequencies(1,ylist)/1000));
    ylabel('Frequency (kHz)');
    cb = colorbar;
    ylabel(cb,'Firing rate (Hz)');
    
    subplot(3,1,3);
    semilogx(frequencies(1,:),tuningCurve);
    xlabel('Frequency (Hz)');
    ylabel('Firing rate (Hz)');
    
    saveas(f1,fullfile('D:\KiloSort\',mouseID,session,folder,'TuningCurves\Figures\',['Unit ' num2str(nData.CellInfo(4)) ' tuning plots.fig']));
    saveas(f1,fullfile('D:\KiloSort\',mouseID,session,folder,'TuningCurves\Figures\',['Unit ' num2str(nData.CellInfo(4)) ' tuning plots.jpg']));
    close(f1);
end

analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.folder = folder;
analysisParams.stimPath = stimPath;
analysisParams.binWindows = binWindows; %ms
analysisParams.analysisWindow = analysisWindow;
save(fullfile('D:\KiloSort\',mouseID,session,folder,'TuningCurves\TuningData.mat'),'unitData','analysisParams');
            