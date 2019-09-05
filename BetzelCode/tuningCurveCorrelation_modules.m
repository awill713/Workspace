clear
% close all
mouse = {'K048','K056','K070','K073','K074'};
% mouse = {'K048','K056'};

for mm = 1:length(mouse)
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    r = find(~isnan(x(:,1))); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    folderNos = x(r,5);
    dataLoc = 'D:\2Pdata\data\';
    
    for ii = 1:length(r)
        disp(ii)
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        
        if ~isempty(fnt)
%              return
            
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
            preEv = proc.preEv;
            fr = exptInfo.fr;
            p_mod = zeros(1,size(proc.raster,3)); bf = p_mod;
            t_resp=zeros(size(proc.raster,3),length(uT));
            for jj = 1:size(proc.raster,3)
                raster = proc.raster(:,:,jj);
                for kk = 1:size(raster,1)
                    raster(kk,:) = SmoothGaus(raster(kk,:),2);
                end
                resp = [];
                for tt = 1:length(uT)
                    rows = find(it==tt);
                    try
                        resp(:,tt) = nanmean(raster(rows,preEv+1:preEv+ceil(fr)),2);
                    catch
                        resp(:,tt) = NaN;
                        resp(1:length(rows),tt) = nanmean(raster(rows,preEv+1:preEv+ceil(fr)),2);
                    end
                end
                p_mod(jj) = kruskalwallis(resp,uT,'off');
                t_resp(jj,:) = nanmean(resp);
                [~,bf(jj)] = max(mean(resp));
            end
            
            cellModules = netanal_tone.conCom;
            tunedCellModules = netanal_tone.conCom(p_mod<0.05,:);
            
           [rho,p_corr] = corr(t_resp');
           [~,index]=sort(cellModules(:,1));
           imagesc(rho(index,index))
           colormap gray
           uMod = unique(cellModules(:,1));
           for kk = 1:length(uMod)
               nMod(kk) = sum(netanal_tone.conCom(:,1)==uMod(kk));
               cs = cumsum(nMod(1:kk));
               if kk==1
                   outline = [0.5,cs(end)+0.5];
               else
                   outline = [cs(end-1)+0.5,cs(end)+0.5];
               end
               patch([outline(1) outline(2) outline(2) outline(1)],[outline(2) outline(2) outline(1) outline(1)],'k','facealpha',0,...
                   'edgeColor','k','LineWidth',1)
           end
           
           for jj = 1:size(cellModules,2)
               cwm = NaN(length(rho),1); com = cwm;
               for kk = 1:length(rho)
                   
                   kmod = cellModules(kk,jj);
                   
                   if sum(cellModules(:,jj)==kmod)>5
                   
                    corWithinMod = rho(cellModules(:,jj)==kmod,kk);
                    corWithinMod(corWithinMod==1)=[];
                    corOutsideMod = rho(cellModules(:,jj)~=kmod,kk);
                   
                    cwm(kk) = mean(corWithinMod);
                    com(kk) = mean(corOutsideMod);
                   end
               end
           end
               
               
           
                
                
                
                
                
                
            
            
            
            
            
            
            
