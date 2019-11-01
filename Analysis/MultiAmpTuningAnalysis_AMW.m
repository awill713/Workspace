
clear;

mouseID = 'AW000';
session = 'Session1';
folder = '2019-06-19_18-56-06';
stimPath = 'D:\Code\TuningCurve\TuningCurveMultiAmpSpkrLeft1-70k';
analysisWindow = [-50 200];
quantWindow = [10 30]; %ms after stimulus onset
rasterHeatMap = 'jet';

dataFolder = fullfile('D:\KiloSort\',mouseID,session,folder,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'TuningCurveMultiAmpSpkr*'));

newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'MultiAmpResponses','Figures');
mkdir(newDir);

load(stimPath);
frequencies = unique(STIM.freqOrder);
amplitudes = unique(STIM.ampOrder);
uniqueTones = length(frequencies);
uniqueAmplitudes = length(amplitudes);

totalUnits = length(dataFiles);

analysisEdges = analysisWindow*1000;
quantificationEdges = quantWindow*1000; %temporal resolution of cheetah acquisition is microseconds
frScalar = 1000/(quantWindow(2)-quantWindow(1));
map = colormap(rasterHeatMap);close;

for n = 1:totalUnits
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    tuningRaster = cell(uniqueTones,uniqueAmplitudes,3); %spike times (1), trial number (y-axis) (2), color (3)
    tuningMat = zeros(uniqueTones,uniqueAmplitudes);
    for f = 1:length(frequencies)
        freq = frequencies(f);
        freqEvents = find(STIM.freqOrder==freq);
        color = map(round(f/length(frequencies)*length(map)),:);
        for a = 1:length(amplitudes)
            amp = amplitudes(a);
            ampEvents = find(STIM.ampOrder==amp);
            event1 = intersect(freqEvents,ampEvents);
            
            eventsOfInterest = event1:(uniqueTones*uniqueAmplitudes):length(eventTimes);
            
            spikeCount = zeros(length(eventsOfInterest),1);
            rasterX = [];
            rasterY = [];
            tempY = (a-1)*length(eventTimes)/uniqueAmplitudes + (f-1)*length(eventTimes)/(uniqueAmplitudes*uniqueTones) +1;
            for e = 1:length(eventsOfInterest)
                time = eventTimes(eventsOfInterest(e));
                quant = length(find(spikeTimes>(time+quantificationEdges(1)) & spikeTimes<(time+quantificationEdges(2))));
                spikeCount(e,1) = quant;
                
                tempY = tempY+1;
                alignedSpikes = spikeTimes-time;
                trialSpikes = alignedSpikes(find(alignedSpikes>analysisEdges(1) & alignedSpikes<analysisEdges(end)))/1000;
                rasterX = [rasterX trialSpikes];
                rasterY = [rasterY repmat(tempY,[1 length(trialSpikes)])];
            end
        
        tuningRaster{f,a,1} = rasterX;
        tuningRaster{f,a,2} = rasterY;
        tuningRaster{f,a,3} = color;
        tuningMat(f,a) = mean(spikeCount)*frScalar;
        end
    end
    
    unitData(n).raster = tuningRaster;
    unitData(n).freqAmpTuning = tuningMat;
    unitData(n).type = nData.CellInfo(6);
    
    f1 = figure;
    subplot(2,1,1); colormap 'jet';
    for r = 1:size(tuningRaster,1)
        for col = 1:size(tuningRaster,2)
            scatter(tuningRaster{r,col,1},tuningRaster{r,col,2},[],tuningRaster{r,col,3},'.');
            hold on;
        end
    end
    for a = 1:uniqueAmplitudes-1
        line([analysisWindow(1) analysisWindow(end)],[a*length(eventTimes)/uniqueAmplitudes a*length(eventTimes)/uniqueAmplitudes],'LineStyle',':');
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial number and type');
    c = colorbar;
    set(c,'YTick',[0 1]);
    set(c,'YTickLabel',[frequencies(1,1) frequencies(1,end)])
    ylabel(c,'Frequency');
    title(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6))]);
    
    subplot(2,1,2);
    imagesc(flip(tuningMat'));
    xlist = [1 10:10:50];
    xticks(xlist);
    xticklabels(frequencies(xlist)/1000)
    xlabel('Frequency (kHz)');
    ylist = [1:2:uniqueAmplitudes];
    yticks(ylist);
    yticklabels(flip(20*log10(amplitudes(ylist)*10)+70));
    ylabel('SPL (dB)');
    cb = colorbar;
    ylabel(cb,'Firing rate (Hz)');

    
    saveas(f1,fullfile(newDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.fig']));
    saveas(f1,fullfile(newDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.jpg']));
    close(f1);
end

analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.folder = folder;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow; %ms
analysisParams.quantificationWindow = quantWindow;
save(fullfile('D:\KiloSort\',mouseID,session,folder,'MultiAmpResponses\TuningData.mat'),'unitData','analysisParams');
            