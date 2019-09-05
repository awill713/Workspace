clear
% close all
mouse = {'K070','K073','K074'};
sparseCells = load('D:\Github\Kath\Kath_2Panalysis\fear\sparseGoodCells.mat');
figure('position',[1          41        1920         963])
for mm = 3:length(mouse)
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    r = find(x(:,2)==1); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    folderNos = x(r,5);
    dataLoc = 'D:\2Pdata\data\';
    
    for ii = 1:length(r)
        disp(ii)
        
        fns = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        %             load([fnt.folder filesep fnt.name])
        load([fns.folder filesep fns.name])
        
        for jj = 1:size(netanal_tone.conCom,2)
            hold on
            subplot(2,3,1)
            uMod = unique(netanal_tone.conCom(:,jj));
            nMod = zeros(1,length(uMod));
            cm = colormap('hsv');
            if length(cm)<length(uMod)
                cm = [cm;colormap('summer')];
            end
            cm = cm(1:floor(length(cm)/length(uMod)):floor(length(cm)/length(uMod))*length(uMod),:);
            
            index = fcn_order_partition(netanal_tone.rho,netanal_tone.conCom(:,jj));
            imagesc(netanal_tone.rho(index,index))
            colormap gray
            title(['# of modules: ' num2str(length(uMod))]);
            hold on
            
            for kk = 1:length(uMod)
                nMod(kk) = sum(netanal_tone.conCom(:,jj)==uMod(kk));
                cs = cumsum(nMod(1:kk));
                if kk==1
                    outline = [0.5,cs(end)+0.5];
                else
                    outline = [cs(end-1)+0.5,cs(end)+0.5];
                end
                patch([outline(1) outline(2) outline(2) outline(1)],[outline(2) outline(2) outline(1) outline(1)],cm(kk,:),'facealpha',0,...
                    'edgeColor',cm(kk,:),'LineWidth',1)
            end
            axis equal tight
            
            subplot(2,3,4)
            imagesc(brighten(spatialInfo.im,1))
            colormap gray
            axis equal tight off
            hold on
            
            for kk = 1:size(netanal_tone.conCom,1)
                patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),cm(netanal_tone.conCom(kk,jj),:),'edgealpha',0,'facealpha',0.5)
            end
            
            %%
            mod2plot = find(nMod>10);
            indd = 0;
            for kk = 1:length(mod2plot)
                c2p = find(netanal_tone.conCom(:,jj)==mod2plot(kk));
                c2p = c2p(stats.driven(c2p)<0.05);
                indd(kk+1) = length(c2p)/10;
                spInd = intersect(find(contains(sparseCells.d(:,1),exptInfo.mouse)),intersect([sparseCells.d{:,2}],spatialInfo.masterCellIndex(c2p,2)));
                [uT,~,ict] = unique(proc.stimOrder(:,1));
                fresp = zeros(length(c2p),length(uT));
                fstd = fresp;
                for nn = 1:length(c2p)
                    for tt = 1:length(uT)
                        fresp(nn,tt) = mean(mean(proc.raster(ict==tt,proc.preEv+1:proc.preEv+floor(exptInfo.fr),c2p(nn)),2));
                        %                              fstd(nn,tt) = std(mean(proc.raster(ict==tt,proc.preEv+1:proc.preEv+floor(exptInfo.fr),c2p(nn)),2));
                    end
                end
                subplot(2,3,[2 5])
                [max_fresp,maxind] = max(fresp,[],2);
                [min_fresp] = min(fresp,[],2);
                xw = fresp-repmat(min_fresp,1,size(fresp,2));
                x = xw./repmat(max(xw,[],2),1,size(fresp,2));
                for nn = 1:size(x,1)
                    x(nn,:) = SmoothGaus(x(nn,:),1);
                end
                [~,sort_x] = sort(maxind);
                a = cumsum(indd(1:kk));
                m2pind = a(end)+0.1:0.1:a(end)+length(c2p)/10;
                plot((x(sort_x,:)+(m2pind)')','-','LineWidth',1,'color',cm(mod2plot(kk),:))
                hold on
                set(gca,'Ytick',0:18:m2pind,'YTickLabel',(0:18:m2pind));
                axis tight
                set(gca,'FontSize',14,'XTick',1:length(uT),'XTickLabel',uT)
                ylabel('\DeltaF/F (stds)')
                xlabel('Frequency (Hz)')
                xtickangle(45)
                
                
                subplot(2,3,3)
                plot(mean(x),'-','LineWidth',2,'color',cm(mod2plot(kk),:))
                hold on
                patch([1:length(uT) length(uT):-1:1],[mean(x,1)+(std(x,[],1)/sqrt(size(x,1)-1)) fliplr(mean(x,1)-(std(x,[],1)/sqrt(size(x,1)-1)))],cm(mod2plot(kk),:),'FaceAlpha',0.3,'edgealpha',0)
                plot([1 length(uT)],[0 0],'k-')
            end
            axis tight
            set(gca,'FontSize',14,'XTick',1:length(uT),'XTickLabel',uT)
            ylabel('\DeltaF/F (stds)')
            xlabel('Frequency (Hz)')
            xtickangle(45)
            title([mouse(mm) ' ' fns.name])
            
            %%
%             pause()
%             clf
%             print(gcf, '-depsc', '-r400', ['D:\2Pdata\networkFigs\' fns.name(1:end-4) '_meanFreqResp.eps'])
            saveFigPDF(gcf,[ 1920         963],['D:\2Pdata\networkFigs\' fns.name(1:end-4) '_meanFreqResp.pdf'])
%             print(gcf, '-dpdf', '-r400','-bestfit', ['D:\2Pdata\networkFigs\' fns.name(1:end-4) '_meanFreqResp.pdf'])
        clf
        end
                  
    end
end
