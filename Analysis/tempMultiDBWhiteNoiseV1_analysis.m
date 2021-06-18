

dataPaths{1} = fullfile('EP013','AW182','20210428');
dataPaths{2} = fullfile('EP013','AW183','20210428');
dataPaths{3} = fullfile('EP013','AW184','20210429');


intensity = [0 10 30 50 70];

high = zeros(0,5);
highPSTH = zeros(0,1291);
highNoStim = zeros(0,1291);
low = zeros(0,5);
lowPSTH = zeros(0,1291);
lowNoStim = zeros(0,1291);

for dp = 1: length(dataPaths)
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'WhiteNoiseMultidB','WhiteNoiseMultidB_DATA.mat');
    load(dataFile);
    
    
    for u =1:length(unitData)
        
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)
            
            if unitData(u).meanResponse(5,4) > unitData(u).meanResponse(1,4)
                high = [high; unitData(u).meanResponse(:,4)'];
                highPSTH = [highPSTH; unitData(u).frTrain(5,:)];
                highNoStim = [highNoStim; unitData(u).frTrain(1,:)];
            else
                low = [low; unitData(u).meanResponse(:,4)'];
                lowPSTH = [lowPSTH; unitData(u).frTrain(5,:)];
                lowNoStim = [lowNoStim; unitData(u).frTrain(1,:)];
            end
        end
    end
end
     
%%
binEdges = analysisParams.analysisWindow(1)+analysisParams.frBinWidth:1:analysisParams.analysisWindow(2);


figure;
hCurve = mean(high); hCurveSTD = std(high)/sqrt(size(high,1));
hP = mean(highPSTH); hPSTD = std(highPSTH)/sqrt(size(highPSTH,1));
hPN = mean(highNoStim); hPNSTD = std(highNoStim)/sqrt(size(highNoStim,1));

lCurve = mean(low); lCurveSTD = std(low)/sqrt(size(low,1));
lP = mean(lowPSTH); lPSTD = std(lowPSTH)/sqrt(size(lowPSTH,1));
lPN = mean(lowNoStim); lPNSTD = std(lowNoStim)/sqrt(size(lowNoStim,1));


subplot(2,2,1);
hold on;
hA = area(intensity,[hCurve-hCurveSTD; 2*hCurveSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(intensity,hCurve,'Color',[0 0 0]);
xlabel('Intensity (dB)');
ylabel('Firing rate (Hz)');

subplot(2,2,2);
hold on;
hA = area(binEdges,[hPN-hPNSTD; 2*hPNSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0.5 0.5 0.5];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges,[hP-hPSTD; 2*hPSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges,hPN,'Color',[0.5 0.5 0.5]);
plot(binEdges,hP,'Color',[0 0 0]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');

subplot(2,2,3);
hold on;
hA = area(intensity,[lCurve-lCurveSTD; 2*lCurveSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(intensity,lCurve,'Color',[0 0 0]);
xlabel('Intensity (dB)');
ylabel('Firing rate (Hz)');

subplot(2,2,4);
hold on;
hA = area(binEdges,[lPN-lPNSTD; 2*lPNSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0.5 0.5 0.5];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges,[lP-lPSTD; 2*lPSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges,lPN,'Color',[0.5 0.5 0.5]);
plot(binEdges,lP,'Color',[0 0 0]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');