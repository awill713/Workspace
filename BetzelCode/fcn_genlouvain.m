function [ci,q] = fcn_genlouvain(R,gammavals,nreps)
N = length(R);
ngam = length(gammavals);
ci = zeros(N,nreps,ngam);
q = zeros(nreps,ngam);
for igam = 1:ngam
    fprintf('gamma %i/%i',igam,ngam);
    tic;
    B = (R - gammavals(igam)).*~eye(N);
    for irep = 1:nreps
        [ci(:,irep,igam),q(irep,igam)] = genlouvain(B);
    end
    fprintf(' ... %.2f s\n',toc);
end