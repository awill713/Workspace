function ciu = fcn_consensus_communities(ci,niter,vis)
if ~exist('vis','var')
    vis = false;
end
N = size(ci,1);
mask = triu(ones(N),1) > 0;
d = agreement(ci);
goFlag = length(unique(d));
if goFlag <= 2
    CiCon = ci;
end
while goFlag > 2
    
    mu = mean(d(mask));
    b = d - mu;
    b(1:(N + 1):end) = 0;
    CiCon = zeros(N,niter);
    for iRep = 1:niter
        CiCon(:,iRep) = genlouvain(b);
        if vis
            imagesc(CiCon); drawnow;
        end
    end
    d = agreement(CiCon);
    goFlag = length(unique(d));
    
end
ciu = CiCon(:,1);