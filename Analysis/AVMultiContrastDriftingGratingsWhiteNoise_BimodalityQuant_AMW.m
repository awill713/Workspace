
clear all;
% close all;

% saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
saveDir = fullfile('D:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','D prime - individual neuron');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP004','AW117','20200221-1');
dataPaths{2} = fullfile('EP004','AW117','20200221-2');
dataPaths{3} = fullfile('EP004','AW118','20200221-1');
dataPaths{4} = fullfile('EP004','AW118','20200221-2');
dataPaths{5} = fullfile('EP004','AW121','20200226-1');
dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{7} = fullfile('EP004','AW124','20200303-1');
dataPaths{8} = fullfile('EP004','AW124','20200303-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


count = 0;
for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Recat','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    
    if ~exist('rVector','var')
        rVector = cell(1,length(contrasts));
        osiMat = cell(1,length(contrasts)); %[osi1V osi2V osi1AV osi2AV dp u neuronNumber]
        dsiMat = cell(1,length(contrasts)); %[dsi1V dsi2V dsi1AV dsi2AV dp u neuronNumber]
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
%         if (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
%                 && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
%                 && (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
%                 && unitData(u).type~=0
            
%                 && (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
        if ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)
            count = count+1;
                            
            responses = unitData(u).meanResponse(:,4);
            
            for c = 1:length(contrasts)
                
                vResp = responses((c-1)*length(orientations)+1:c*length(orientations));
                avResp = responses((c-1)*length(orientations)+1+length(orientations)*length(contrasts):c*length(orientations)+length(contrasts)*length(orientations));
                
                vRCos = sum(vResp'.*cosd(orientations)) / sum(vResp);
                vRSin = sum(vResp'.*sind(orientations)) / sum(vResp);
                avRCos = sum(avResp'.*cosd(orientations)) / sum(avResp);
                avRSin = sum(avResp'.*sind(orientations)) / sum(avResp);
                
                vR = sqrt(vRCos^2 + vRSin^2);
                avR = sqrt(avRCos^2 + avRSin^2);
                rVector{c} = [rVector{c}; vR avR dp u neuronNumber];
                
                [vPref, vPrefInd] = max(vResp);
                [avPref, avPrefInd] = max(avResp);
                
                orthoIndV = mod([vPrefInd+length(orientations)/4 vPrefInd-length(orientations)/4]-1,length(orientations))+1;
                orthoIndAV = mod([avPrefInd+length(orientations)/4 avPrefInd-length(orientations)/4]-1,length(orientations))+1;
                oppIndV = mod(vPrefInd+length(orientations)/2-1,length(orientations))+1;
                oppIndAV = mod(avPrefInd+length(orientations)/2-1,length(orientations))+1;
                
                osi1V = 1 - mean(vResp(orthoIndV))/vPref;
                osi2V = (vPref - mean(vResp(orthoIndV))) / (vPref + mean(vResp(orthoIndV)));
                osi1AV = 1 - mean(avResp(orthoIndAV))/avPref;
                osi2AV = (avPref - mean(avResp(orthoIndAV))) / (avPref + mean(avResp(orthoIndAV)));
                
                osiMat{c} = [osiMat{c}; osi1V osi2V osi1AV osi2AV dp u neuronNumber];
                
                dsi1V = 1 - vResp(oppIndV)/vPref;
                dsi2V = (vPref - vResp(oppIndV)) / (vPref + vResp(oppIndV));
                dsi1AV = 1 - avResp(oppIndAV)/avPref;
                dsi2AV = (avPref - avResp(oppIndAV)) / (avPref + avResp(oppIndAV));
                
                dsiMat{c} = [dsiMat{c}; dsi1V dsi2V dsi1AV dsi2AV dp u neuronNumber];
            end
        end
    end
end

count

figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    histogram(rVector{c}(:,1),0:0.05:0.6);
    histogram(rVector{c}(:,2),0:0.05:0.6);
    osiDiff = rVector{c}(:,2) - rVector{c}(:,1);
    title(['avg = ' num2str(mean(osiDiff))]);
end
suptitle('rVector');

% figure;
% histogram(rVector{5}(:,1));
% xlabel('||r||');

figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    histogram(osiMat{c}(:,1));
    histogram(osiMat{c}(:,3));
    dsiDiff = osiMat{c}(:,3) - osiMat{c}(:,1);
    title(['avg = ' num2str(mean(dsiDiff))]);
end
suptitle('OSI');

figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    histogram(dsiMat{c}(:,1));
    histogram(dsiMat{c}(:,3));
    dsiDiff = dsiMat{c}(:,3) - dsiMat{c}(:,1);
    title(['avg = ' num2str(mean(dsiDiff))]);
end
suptitle('DSI');

% figure;
% for c = 1:length(contrasts)
%     subplot(1,length(contrasts),c);
%     scatter(dsiMat{c}(:,1),osiMat{c}(:,1));
% end

% figure;
% histogram(rVector{5}(:,1));
% xlabel('||r||');


dPrimeFolder = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','D prime - individual neuron');
dPrimeData = load(fullfile(dPrimeFolder,'dPrimeData.mat'));

figure;
for c = 1:length(contrasts)
    orthoPrimes = dPrimeData.orthoDPrime{c};
    osis = osiMat{c};    
    
    xOSI = osis(ismember(osis(:,5:7),orthoPrimes(:,4:6),'rows'),1);
    yDPrimes = orthoPrimes(:,3);
    
    vectorGroups = discretize(xOSI,0:0.1:1);
    meanDPrime = grpstats(yDPrimes,vectorGroups,'mean')';
    meanOSI = grpstats(xOSI,vectorGroups,'mean')';
    semDPrime = grpstats(yDPrimes,vectorGroups,'sem')';
    cc = grpstats(yDPrimes,vectorGroups,'numel');
    
    subplot(1,length(contrasts),c);hold on;
    scatter(xOSI,yDPrimes);
    
    
    hA = area(meanOSI,[meanDPrime-semDPrime; 2*semDPrime]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(meanOSI,meanDPrime,'Color',[0 0 0]);
    %     ylim([-0.3 0.3]);
    %     xticks(contrasts);
    xlabel('OSI');
    ylabel('\Delta D prime ortho');
    title(['C = ' num2str(c)]);

end


figure;
for c = 1:length(contrasts)
    oppPrimes = dPrimeData.oppDPrime{c};
    dsis = dsiMat{c};    
    
    xDSI = dsis(ismember(dsis(:,5:7),oppPrimes(:,4:6),'rows'),1);
    yDPrimes = oppPrimes(:,3);
    
    vectorGroups = discretize(xDSI,0:0.1:1);
    meanDPrime = grpstats(yDPrimes,vectorGroups,'mean')';
    meanDSI = grpstats(xDSI,vectorGroups,'mean')';
    semDPrime = grpstats(yDPrimes,vectorGroups,'sem')';
    cc = grpstats(yDPrimes,vectorGroups,'numel');
    
    subplot(1,length(contrasts),c);hold on;
    scatter(xDSI,yDPrimes);
    
    
    hA = area(meanDSI,[meanDPrime-semDPrime; 2*semDPrime]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(meanDSI,meanDPrime,'Color',[0 0 0]);
    %     ylim([-0.3 0.3]);
    %     xticks(contrasts);
    xlabel('DSI');
    ylabel('\Delta D prime opposite');
    title(['C = ' num2str(c)]);

end

% save(fullfile(saveDir,'SelectivityIndices.mat'),'osiMat','dsiMat');