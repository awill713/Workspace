clear

load('D:\2Pdata\communityData.mat')
mouse = {'K056','K070','K073','K074'};


for mm = 1:length(mouse)
    
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
    
    for ii = 10:length(c) % for each recording...
        
        disp(ii)
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        load([fnt.folder filesep fnt.name])
        
        if iscell(stimInfo)
            stimInfo = stimInfo{1};
        end
        order = stimInfo.index(stimInfo.order,:);
        if length(order)<size(proc.raster,1)
            order = repmat(order,ceil(size(proc.raster,1)/length(order)),1);
            order = order(1:size(proc.raster,1),:);
        end
        
        [uT,~,it] = unique(order(:,1));
        [uA,~,ia] = unique(order(:,2));
        preEv = proc.preEv;
        fr = exptInfo.fr;
        p_mod = zeros(1,size(proc.raster,3)); bf = p_mod; sig = p_mod;
        mresp = zeros(size(proc.raster,3),length(uT));
        % is the cell driven?
        for jj = 1:size(proc.raster,3)
            raster = proc.raster(:,:,jj);
            mca = zeros(length(uA),length(uT));
            h = zeros(length(uA),length(uT)); p = h;
            for aa = 1:length(uA)
                for tt = 1:length(uT)
                    rows = find(it==tt & ia==aa);
                    mca(aa,tt) = mean(nanmean(raster(rows,preEv+1:preEv+ceil(fr))));
                    [h(aa,tt),p(aa,tt)] = ttest(mean(raster(rows,1:preEv),2),mean(raster(rows,preEv+1:preEv+ceil(fr)),2));
                end
            end
            
            [h,p] = ttest(mean(raster(:,1:preEv),2),mean(raster(:,preEv+1:preEv+ceil(fr)),2));
            
            mresp(jj,:) = mean(mca,1);
            pmod = kruskalwallis(nanmean(raster(:,preEv+1:preEv+ceil(fr)),2),it,'off');
            if pmod<0.05
                %             if p<0.05
                sig(jj) = 1;
            end
        end
        
        comInd = find(strcmp(data(mm).days,num2str(dates(ii))));
        drivenCellModules = data(mm).communities{comInd}(sig==1,:); % which modules are the driven neurons in?
        
        % How many of ANY cells are in each module? And how many driven
        % in each module?
        uModules = cell(1,size(data(mm).communities{comInd},2));
        n = uModules;
        n_driven = uModules;
        for jj = 1:size(data(mm).communities{comInd},2)
            uModules{jj} = unique(data(mm).communities{comInd}(:,jj)); % what are the unique modules?
            for kk = 1:length(uModules{jj})
                n{jj}(kk) = sum(data(mm).communities{comInd}(:,jj)==uModules{jj}(kk));
                n_driven{jj}(kk) = sum(data(mm).communities{comInd}(sig==1,jj)==uModules{jj}(kk));
            end
        end
        
        micronsPerPixel = 1.26890590123792;
        sigCells = find(sig==1);
        for jj = 13%:length(n)
            %             figure
%             imagesc(spatialInfo.im)
            hold on
            c = colormap('jet');
            c = c(1:floor(size(c,1)/length(uModules{jj})):floor(size(c,1)/length(uModules{jj}))*length(uModules{jj}),:);
            c2 = colormap('jet');
            c2 = c(1:floor(size(c,1)/sum(sig)):floor(size(c,1)/sum(sig))*sum(sig),:);
            colormap gray
            for kk = 1:size(proc.raster,3)
                patch(spatialInfo.ROIs{kk}(:,1)/micronsPerPixel,spatialInfo.ROIs{kk}(:,2)/micronsPerPixel,c(data(mm).communities{comInd}(kk,jj),:),'facealpha',0.5,'edgealpha',1);
                if ismember(kk,sigCells)
                    plot(spatialInfo.ROIs{kk}(:,1)/micronsPerPixel,spatialInfo.ROIs{kk}(:,2)/micronsPerPixel,'k','LineWidth',2)
                    [~, bf] = max(mresp(kk,:));
                    bf = uT(bf);
                    text(spatialInfo.centroid(kk,2)/micronsPerPixel,spatialInfo.centroid(kk,1)/micronsPerPixel+12/micronsPerPixel,num2str(round(bf/1000,2,'significant')),'Color','k','FontSize',8)
                end
            end
            
            pause()
            clf
        end

        axis tight equal
        set(gca,'TickDir','out','FontSize',14)
        ylabel('Dorsal-ventral (\mum)')
        xlabel('Rostral-caudal (\mum)')

        
        micronsPerPixel = 1.26890590123792;
        sigCells = find(sig==1);
        for jj = 13%:length(n)
            %             figure
            %             imagesc(spatialInfo.im)
            hold on
            c = colormap('jet');
            c = c(1:floor(size(c,1)/length(uModules{jj})):floor(size(c,1)/length(uModules{jj}))*length(uModules{jj}),:);
            c2 = colormap('jet');
            c2 = c(1:floor(size(c,1)/sum(sig)):floor(size(c,1)/sum(sig))*sum(sig),:);
            colormap gray
            for kk = 1:size(proc.raster,3)
                hold on
                t = linspace(0, 2*pi);
                r = 6/micronsPerPixel;
                x = r*cos(t)+spatialInfo.centroid(kk,1)/micronsPerPixel;
                y = r*sin(t)+spatialInfo.centroid(kk,2)/micronsPerPixel;
                figure(1)
                patch(y, x, c(data(mm).communities{comInd}(kk,jj),:),'facealpha',0.5,'edgealpha',0);
                hold on               
            end
            
%             pause()
%             clf
        end

        axis tight equal
        set(gca,'TickDir','out','FontSize',14)
        ylabel('Dorsal-ventral (\mum)')
        xlabel('Rostral-caudal (\mum)')
        
    end
end

