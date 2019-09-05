clear
% close all
mouse = {'K056','K070','K073','K074'};
load('D:\2Pdata\communityData.mat')
% figure('position',[1          41        1920         963])

for mm = 2%:length(mouse)
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    recs = x(:,3);
    r = ismember(recs,str2double(data(mm).days)); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    folderNos = x(r,5);
    dataLoc = 'D:\2Pdata\data\';
    
    if strcmp(mouse{mm},'K070')
        folderNos(ismember(dates,20170731))=[];
        dates(ismember(dates,20170731))=[];
        folderNos(ismember(dates,20170812))=[];
        dates(ismember(dates,20170812))=[];
    end
    f1 = figure;
    f2 = figure;
    for ii = 10 %1:length(dates) % for each recording...
        
        disp(ii)
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        load([fnt.folder filesep fnt.name])
        
        comInd = find(strcmp(data(mm).days,num2str(dates(ii))));
        
        % How many of ANY cells are in each module? And how many driven
        % in each module?
        gamma = 13;
        modules = data(mm).communities{comInd}(:,gamma);
        uMod = unique(modules);
        for jj = 1:length(uMod)
            nMod(jj) = sum(modules==uMod(jj));
        end
        cm = colormap('jet');
        cm = cm(1:floor(length(cm)/length(uMod)):floor(length(cm)/length(uMod))*length(uMod),:);
        mod2plot = find(nMod>3);
        indd = 0;
        for kk = 1:length(mod2plot)
            c2p = find(modules==mod2plot(kk));
            c2p = c2p(stats.driven(c2p)<0.05);
            indd(kk+1) = length(c2p)/10;
%             spInd = intersect(find(contains(sparseCells.d(:,1),exptInfo.mouse)),intersect([sparseCells.d{:,2}],spatialInfo.masterCellIndex(c2p,2)));
            [uT,~,ict] = unique(proc.stimOrder(:,1));
            fresp = zeros(length(c2p),length(uT));
            fstd = fresp;
            for nn = 1:length(c2p)
                for tt = 1:length(uT)
                    fresp(nn,tt) = mean(mean(proc.raster(ict==tt,proc.preEv+1:proc.preEv+floor(exptInfo.fr),c2p(nn)),2));
                    %                              fstd(nn,tt) = std(mean(proc.raster(ict==tt,proc.preEv+1:proc.preEv+floor(exptInfo.fr),c2p(nn)),2));
                end
            end
            figure(f1)
            [max_fresp,maxind] = max(fresp,[],2);
            [min_fresp] = min(fresp,[],2);
            xw = fresp-repmat(min_fresp,1,size(fresp,2));
            x = xw./repmat(max(xw,[],2),1,size(fresp,2));
            for nn = 1:size(x,1)
                x(nn,:) = SmoothGaus(x(nn,:),1);
                x(nn,:) = (x(nn,:)-min(x(nn,:)))/max((x(nn,:)-min(x(nn,:))));
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
            
            
            figure(f2)
            plot(mean(x),'-','LineWidth',2,'color',cm(mod2plot(kk),:))
            hold on
            errorbar(mean(x,1),std(x,[],1)/sqrt(size(x,1)-1),'LineWidth',2,'LineStyle','none','Color',cm(mod2plot(kk),:))
%             patch([1:length(uT) length(uT):-1:1],[mean(x,1)+(std(x,[],1)/sqrt(size(x,1)-1)) fliplr(mean(x,1)-(std(x,[],1)/sqrt(size(x,1)-1)))],cm(mod2plot(kk),:),'FaceAlpha',0.3,'edgealpha',0)
            plot([1 length(uT)],[0 0],'k-')
            axis tight
        set(gca,'FontSize',14,'XTick',1:length(uT),'XTickLabel',uT)
        ylabel('\DeltaF/F (stds)')
        xlabel('Frequency (Hz)')
        xtickangle(45)
        title([mouse(mm) ' ' fnt.name])
        end

    end
    
end
% end
