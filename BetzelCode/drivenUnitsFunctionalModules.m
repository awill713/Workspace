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
            
            fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
            load([fnt.folder filesep fnt.name])
            load([fns.folder filesep fns.name],'netanal')
            
            if ~isfield(netanal,'shuf_sig_driven')
                
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
                for jj = 1:size(proc.raster,3)
                    raster = proc.raster(:,:,jj);
                    for kk = 1:size(raster,1)
                        raster(kk,:) = SmoothGaus(raster(kk,:),2);
                    end
                    
                    mca = zeros(length(uA),length(uT));
                    h = zeros(length(uA),length(uT)); p = h;
                    for aa = 1:length(uA)
                        for tt = 1:length(uT)
                            rows = find(it==tt & ia==aa);
                            mca(aa,tt) = mean(nanmean(raster(rows,preEv+1:preEv+ceil(fr))));
                            [h(aa,tt),p(aa,tt)] = ttest(mean(raster(rows,1:preEv),2),mean(raster(rows,preEv+1:preEv+ceil(fr)),2));
                        end
                    end
                    
                    
                    if sum(h(:))>2
                        sig(jj)=1;
                    end
%                     p_mod(jj) = kruskalwallis(resp,uT,'off');
%                     [~,bf(jj)] = max(mean(resp));
                end
                
                drivenCellModules = netanal.conCom(sig==1,:); % which modules are the tuned neurons in?
                
                % How many of ANY cells are in each module? And how many tuned
                % in each module?
                uModules = cell(1,size(netanal.conCom,2));
                n = uModules;
                n_tuned = uModules;
                for jj = 1:size(netanal.conCom,2)
                    uModules{jj} = unique(netanal.conCom(:,jj)); % what are the unique modules?
                    for kk = 1:length(uModules{jj})
                        n{jj}(kk) = sum(netanal.conCom(:,jj)==uModules{jj}(kk));
                        n_tuned{jj}(kk) = sum(netanal.conCom(sig==1,jj)==uModules{jj}(kk));
                    end
                end
                
                % Now perform a shuffle and randomly select the number of tuned
                % cells from the modules you have and see if any are
                % significantly different than chance (2.5/97.5%)
                nShuf=2000;
                n_mod_shuf = cell(1,length(uModules));
                n_shuf_sig = n_mod_shuf;
                
                for jj = 1:length(uModules)
                    n_shuf = zeros(nShuf,length(drivenCellModules));
                    for ss = 1:nShuf
                        s=randperm(sum(n{jj}),length(drivenCellModules));
                        for kk = 1:length(s)
                            n_shuf(ss,kk) = find(s(kk)<=cumsum(n{jj}),1,'first');
                        end
                    end
                    
                    for kk = 1:length(uModules{jj})
                        n_mod_shuf{jj}(:,kk) = sum(n_shuf==uModules{jj}(kk),2);
                        n_shuf_sig{jj}(kk) = n_tuned{jj}(kk)<prctile(n_mod_shuf{jj}(:,kk),2.5) | n_tuned{jj}(kk)>prctile(n_mod_shuf{jj}(:,kk),97.5);
                    end
                end
                
                netanal.shuf_sig_driven = n_shuf_sig;
                save([fns.folder filesep fns.name],'netanal','-append')
                save([fnt.folder filesep fnt.name],'netanal','-append')
                for jj = 1:length(n_shuf_sig)
                    disp(n_shuf_sig{jj})
                end
                clear netanal
            end
        end
    end
end
