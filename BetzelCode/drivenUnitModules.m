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
    
    for ii = 1:length(c) % for each recording...
        
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
            
            mresp(jj,:) = mean(mca,1);
            if sum(h(:))>2
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
        
        % Now perform a shuffle and randomly select the number of tuned
        % cells from the modules you have and see if any are
        % significantly different than chance (2.5/97.5%)
        nShuf = 2000;
        n_driv_shuf = cell(1,length(uModules));
        n_shuf_sig = n_driv_shuf;
        
        for jj = 1:length(uModules)
            n_shuf = zeros(nShuf,length(drivenCellModules));
            for ss = 1:nShuf
                s = randperm(sum(n{jj}),length(drivenCellModules));
                for kk = 1:length(s)
                    n_shuf(ss,kk) = find(s(kk) <= cumsum(n{jj}),1,'first');
                end
            end
            
            for kk = 1:length(uModules{jj})
                n_driv_shuf{jj}(:,kk) = sum(n_shuf == uModules{jj}(kk),2);
                n_shuf_sig{jj}(kk) = n_driven{jj}(kk) < prctile(n_driv_shuf{jj}(:,kk),2.5) | n_driven{jj}(kk) > prctile(n_driv_shuf{jj}(:,kk),97.5);
            end
        end
        
        newNet.n_shuf_driv = n_shuf_sig;
        
        save([fnt.folder filesep fnt.name],'newNet','-append')
%         dc  = find(sig==1);
%         for jj = 1:length(n_shuf_sig)
%             s = find(n_shuf_sig{jj}==1);
%             for ss = 1:length(s)
%                 cells = find(drivenCellModules(:,jj) == s(ss));
%                 disp(cells)
%                 for cc = 1:length(cells)
%                     plot(mean(mresp(dc(cells(cc)),:),1))
%                     hold on
%                 end
%                 pause()
%                 clf
%             end
%         end

    end   
end

