clear

load('D:\2Pdata\communityData.mat')
mouse = {'K070','K073','K074'}; % can't do K056 because only have one sessions of spontaneous for that mouse - have many tone recordings
data = data(2:end);

for mm = 2%:length(mouse)
    
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    recs = x(:,3);
    r = ismember(recs,str2double(data(mm).days)) & x(:,1)==0; % dates to analyse
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
    
    bfs=[];
    cents = [];
    bestFs = [];
    
    for ii = 1:length(dates) % for each recording...
        
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
        
        [~, bf] = max(mresp,[],2);    
        sigCells = find(sig==1);
  
%         figure
%         c = colormap('jet');
%         close(gcf)
%         c = c(1:floor(size(c,1)/length(uT)):floor(size(c,1)/length(uT))*length(uT),:);
%         colormap gray
%         imagesc(spatialInfo.im)
%         for jj = 1:length(sigCells)
%             hold on
%             t = linspace(0, 2*pi);
%             r = 50;
%             x = r*cos(t)+spatialInfo.centroid(sigCells(jj),1);
%             y = r*sin(t)+spatialInfo.centroid(sigCells(jj),2);
%             figure(1)
%             patch(x, y, c(bf(jj),:),'facealpha',0.5,'edgealpha',0);
%             hold on            
%         end
%         axis equal tight
        
%         figure; 

        cents = [cents; [spatialInfo.centroid(sigCells,2)]];
        bestFs = [bestFs; uT(bf(sigCells))];
%         plot([spatialInfo.centroid(sigCells,2)],uT(bf(sigCells)),'o')
%         hold on
%         drawnow

        
    end
%     pause()
%     clf
end

plot(cents,bestFs,'o','LineWidth',2)
