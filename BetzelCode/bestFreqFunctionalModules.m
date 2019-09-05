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
        fns = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_spontaneous*.mat']);
        
        if ~isempty(fns)
%             return
            fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
            load([fnt.folder filesep fnt.name])
            load([fns.folder filesep fns.name],'netanal')
            
            if ~isfield(netanal,'sig_bfDistribution')
                
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
                    [~,bf(jj)] = max(mean(resp));
                end
                
                bf = uT(bf);
                
                tunedCellModules = netanal.conCom(p_mod<0.05,:); % which modules are the tuned neurons in?
                tunedBfs = bf(p_mod<0.05);
                uTunedBfs = unique(tunedBfs);
                % How many of ANY cells are in each module? And how many tuned
                % in each module?
                nBfsMod = cell(1,size(netanal.conCom,2));
                uTunedCellModules = nBfsMod;
                for jj = 1:size(netanal.conCom,2)
                    uTunedCellModules{jj} = unique(tunedCellModules(:,jj));
                    for nn = 1:length(uTunedCellModules{jj})
                        for kk = 1:length(uTunedBfs)
                            nBfsMod{jj}(kk,nn) = sum(tunedBfs(tunedCellModules(:,jj)==uTunedCellModules{jj}(nn))==uTunedBfs(kk));
                            %                         n_tuned{jj}(kk) = sum(netanal.conCom(p_mod<0.05,jj)==uModules{jj}(kk));
                        end
                    end
                end
                
                % Now perform a shuffle and randomly select the number of tuned
                % cells from the modules you have and see if any are
                % significantly different than chance (2.5/97.5%)
                nShuf=2000;
                nBfsMod_shuf = cell(1,size(netanal.conCom,2));
                for ss = 1:nShuf
                    bf_shuf = tunedBfs(randperm(length(tunedBfs)));
                    % How many of ANY cells are in each module? And how many tuned
                    % in each module?
                    for jj = 1:size(netanal.conCom,2)
                        uTunedCellModules{jj} = unique(tunedCellModules(:,jj));
                        for nn = 1:length(uTunedCellModules{jj})
                            for kk = 1:length(uTunedBfs)
                                nBfsMod_shuf{jj}(kk,nn,ss) = sum(bf_shuf(tunedCellModules(:,jj)==uTunedCellModules{jj}(nn))==uTunedBfs(kk));
                                %                         n_tuned{jj}(kk) = sum(netanal.conCom(p_mod<0.05,jj)==uModules{jj}(kk));
                            end
                        end
                    end
                end
                
                sig_bfDistribution = cell(1,size(netanal.conCom,2));
                for jj = 1:size(netanal.conCom,2)
                    uTunedCellModules{jj} = unique(tunedCellModules(:,jj));
                    for nn = 1:length(uTunedCellModules{jj})
                        for kk = 1:length(uTunedBfs)
                            
                            sig_bfDistribution{jj}(kk,nn) = nBfsMod{jj}(kk,nn)<prctile(nBfsMod_shuf{jj}(kk,nn,:),2.5)...
                                | nBfsMod{jj}(kk,nn)>prctile(nBfsMod_shuf{jj}(kk,nn,:),97.5);
                            
                        end
                    end
                end
                
                netanal.sig_bfDistribution = sig_bfDistribution;
                save([fns.folder filesep fns.name],'netanal','-append')
%                 save([fnt.folder filesep fnt.name],'netanal','-append')
                for jj = 1:length(sig_bfDistribution)
                    disp(sig_bfDistribution{jj}')
                end
                clear netanal
            end
        end
    end
end
