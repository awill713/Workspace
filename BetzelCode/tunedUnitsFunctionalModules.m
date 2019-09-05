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
        %         fns = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_spontaneous*.mat']);
        
        %         if ~isempty(fns)
        
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        load([fnt.folder filesep fnt.name])
        %             load([fns.folder filesep fns.name],'netanal_tone')
        
        if ~isfield(netanal_tone,'shuf_tuned_sig')
            
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
            
            tunedCellModules = netanal_tone.conCom(p_mod<0.05,:); % which modules are the tuned neurons in?
            
            % How many of ANY cells are in each module? And how many tuned
            % in each module?
            uModules = cell(1,size(netanal_tone.conCom,2));
            n = uModules;
            n_tuned = uModules;
            for jj = 1:size(netanal_tone.conCom,2)
                uModules{jj} = unique(netanal_tone.conCom(:,jj)); % what are the unique modules?
                for kk = 1:length(uModules{jj})
                    n{jj}(kk) = sum(netanal_tone.conCom(:,jj)==uModules{jj}(kk));
                    n_tuned{jj}(kk) = sum(netanal_tone.conCom(p_mod<0.05,jj)==uModules{jj}(kk));
                end
            end
            
            % Now perform a shuffle and randomly select the number of tuned
            % cells from the modules you have and see if any are
            % significantly different than chance (2.5/97.5%)
            nShuf=2000;
            n_mod_shuf = cell(1,length(uModules));
            n_shuf_sig = n_mod_shuf;
            
            for jj = 1:length(uModules)
                n_shuf = zeros(nShuf,length(tunedCellModules));
                for ss = 1:nShuf
                    s=randperm(sum(n{jj}),length(tunedCellModules));
                    for kk = 1:length(s)
                        n_shuf(ss,kk) = find(s(kk)<=cumsum(n{jj}),1,'first');
                    end
                end
                
                for kk = 1:length(uModules{jj})
                    n_mod_shuf{jj}(:,kk) = sum(n_shuf==uModules{jj}(kk),2);
                    n_shuf_sig{jj}(kk) = n_tuned{jj}(kk)<prctile(n_mod_shuf{jj}(:,kk),2.5) | n_tuned{jj}(kk)>prctile(n_mod_shuf{jj}(:,kk),97.5);
                end
            end
            
            netanal_tone.shuf_tuned_sig = n_shuf_sig;
            %                 save([fns.folder filesep fns.name],'netanal_tone','-append')
            save([fnt.folder filesep fnt.name],'netanal_tone','-append')
            for jj = 1:length(n_shuf_sig)
                disp(n_shuf_sig{jj})
            end
            clear netanal_tone
        end
    end
end

