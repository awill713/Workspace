function [dff,mu,sd,gt,z] = fcn_jit(ts,fs,buffsec,nrand)
% FCN_JIT       jitter time series and recompute correlations
%% number of samples/nodes
[T,N] = size(ts);           % number of samples, nodes
if T < N
    ts = ts';
    [T,N] = size(ts);
end
%% estimate correlation
p = round(fs*buffsec);
idx = (p + 1):(T - p);
c = fcn_fisher(corr(ts(idx,:)));
%% initialize some matrices
gt = zeros(N);
cnull = zeros(N,N,nrand);
%% generate nrand jittered ts
for jrand = 1:nrand
    fprintf('randomisation: %d/%d\n',jrand,nrand)
    % +/- some non-zero jitter
    r = randi(p,1,N).*sign(randn(1,N));
    
    % initialize array
    tsr = zeros(T - 2*p,N);
    
    % loop over nodes
    for j = 1:N
        
        % new indices and samples
        jdx = idx + r(j);
        tsr(:,j) = ts(jdx,j);
    end
    
    % compute correlation
    cr = fcn_fisher(corr(tsr));
    
    % was original corr > jittered
    gt = gt + (c > cr);
    
    % store jittered correlation matrix
    cnull(:,:,jrand) = cr;
end
%% post jitter statistics
% mean difference
dff = -mean(bsxfun(@minus,cnull,c),3);

% mean jittered correlation
mu = nanmean(cnull,3);

% standard deviation jittered correlation
sd = nanstd(cnull,[],3);

% zscore correlations
z = (c - mu)./sd;
function rvalsTrans = fcn_fisher(rvals)
%FCN_FISHER         fisher transformation
%
%   RVALSTRANS = FCN_FISHER(RVALS) applies a Fisher transformation to
%                correlation coefficients RVALS.
%
%   Inputs:     RVALS,      raw coefficients
%
%   Outputs:    RVALSTRANS, transformed coefficients
%
%   Richard Betzel, Indiana University, 2012
%

rvalsTrans = 0.5*log((1 + rvals)./(1 - rvals));
rvalsTrans(isinf(rvalsTrans)) = 0;