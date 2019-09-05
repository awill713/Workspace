function [netanal] = networkAssignmentCode_fcn(filename)

% clear all
% close all
addpath(genpath('GenLouvain-master\'))
clc
% filename = 'D:\2Pdata\data\K070_20170817_spontaneous_02.mat';
%% prepare data structure
name = 'dataStruct.mat';
check = dir(name);
if isempty(check)
    load(filename,'calcium','exptInfo'); % change this to load your data
    data.ts = double(calcium.npilSubTraces);
    [COEFF, SCORE] = pca(data.ts');
    data.fr = exptInfo.fr;
    data.proc.gsr = false;
    data.proc.diff = false;
    data.proc.filt.filttype = [];
    data.proc.filt.filtorder = [];
    data.proc.filt.freqrange = [];
    data.net.maxlag_seconds = 5;
%     save(name,'data');
else
    load(name);
end
%% generate network
name = 'exampleNetwork.mat';
check = dir(name);
if isempty(check)
    [R,lag] = fcn_get_network_traces(data);
%   Rho,        neuron x neuron correlation matrix
%   Lag,        neuron x neuron matrix of time lags
%   save(name,'R');
else
    load(name);
end
%% estimate modules
%% modularity maximization
% we've tried a few different clustering approaches. a really popular
% method is "modularity maximization" (newman & girvan 2004).
%
% modularity maximization tries to optimize a particular objective function
% called "modularity" or "Q". the goal is to divide a network into clusters
% whose internal density of connections (mean connection strength among all
% nodes assigned to the same cluster) is maximally greater than a
% chance/null model.
%
% this method has one free parameter, gamma, that determines the
% number/size of clusters. when gamma is small, you get few clusters, but
% when its large you get many smal clusters.
%
% the user also has to specify the chance/random connectivity profile. we
% follow a recent paper and use a "flat" model (bazzi et al 2016, 
% multiscale model.). there are others but this one works well when the 
% network is estimated from a correlation matrix.
%
% the software package we use to implement modularity maximization is the a
% so-called generalized louvain algorithm. you can download it here:
%   http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%
% let me know if you have trouble installing -- i think the newest version
% doesn't require any compiling. just add the directories to your current
% path.
% 
%
% the function "fcn_genlouvain.m" will implement this step.
N = length(R);      % number of nodes
ngam = 41;          % number of gamma values to test
nreps = 100;        % number of restarts of the louvain algorithm
mask = triu(ones(N),1) > 0;     % upper triangle mask
mu = mean(R(mask));             % mean connection weight
sd = std(R(mask));              % stdv connection weight
gammavals = ...                 % gamma values to test
    linspace(0,1,ngam)*sd + mu;

name = 'exampleCommunitiesQmax.mat';
check = dir(name);
if isempty(check)
    ci = zeros(N,nreps,ngam);
    for igam = 1:ngam
        B = (R - gammavals(igam)).*~eye(N); % make a "modularity matrix"
        for irep = 1:nreps
            ci(:,irep,igam) = genlouvain(B);% feed matrix to algorithm
        end
    end
%     save(name,'ci','gammavals');
else 
    load(name);
end
%% choose optimal gamma
% again, lots of strategies for choosing gamma (i.e. number/size of
% clusters). one heuristic is to focus on gamma values for which the
% algorithm repeatedly converges to the same (or similar) cluster
% solutions.
%
% here, we use the z-score of the rand index to measure similarity.
% basically, it tells us how much more similar to cluster partitions are to
% one another compared to a null distribution (traud et al 2011).
%
% we'll focus on local maxima of the median z-score rand index versus gamma
% curve. to identiy ``stable'' maxima, we'll take random sub-samples of the
% all pairwise similarity scores and identify local maxima that are present
% consistently across subsamples.
%
%
% the function "fcn_choose_gamma.m" will implement this step.

subsample = round(nreps/2);             % number of partitions in subsample
nsubsample = 100;                       % number of subsamples
mask = triu(ones(nreps),1) > 0;         % upper triangle mask
masksub = triu(ones(subsample),1) > 0;  % a different upper triangle mask
med = zeros(1,ngam);                    % vector to store median similarity over ALL cluster partitions
medsubsample = zeros(nsubsample,ngam);  % matrix to store median similarity from subsamples
for igam = 1:ngam                       % loop over gamma
    z = zeros(nreps);                   % we'll temporarily store similarity scores here.
    for irep = 1:(nreps - 1)            % loop over pairs of partitions
        for jrep = (irep + 1):nreps
            z(irep,jrep) = fcn_zrand(ci(:,irep,igam),ci(:,jrep,igam));  % compute similarity
            z(jrep,irep) = z(irep,jrep);                                % z-score rand index is symmetric, so z(i,j) is the same as z(j,i)
        end
    end
    med(igam) = median(z(mask));        % get median similarity
    for isub = 1:nsubsample             % randomly subsample
        r = randperm(nreps,subsample);  % select pairs at random
        zsub = z(r,r);                  % extract subsampled matrix
        medsubsample(isub,igam) = median(zsub(masksub));    % compute median z-score rand index
    end
    ph = plot(1:igam,medsubsample(:,1:igam),1:igam,med(1:igam));    % make a plot
    set(ph(end),'color','k','linewidth',2);
    drawnow;
end

pks = zeros(size(medsubsample));
for isub = 1:nsubsample                 % loop over subsamples and find peaks
    idx = findpeaks(medsubsample(isub,:));
    pks(isub,idx.loc) = 1;
end
% mu = mean(pks,1) >= 0.95;               % keep peaks that are consistent in 95% of subsamples
threshold = 0.95;
mu = mean(pks,1) >= threshold;               % keep peaks that are consistent in 95% of subsamples
while sum(mu)==0
    threshold = threshold-0.01;
    disp(threshold)
    mu = mean(pks,1) >= threshold;
end
t = threshold;
while sum(mu)==1
    threshold = threshold-0.01;
    disp(threshold)
    mu = mean(pks,1) >= threshold;
end
t(2) = threshold;
cipks = ci(:,:,mu);
cicon = zeros(N,size(cipks,3));         % do consensus clustering to get "representative" clusters (lancichinetti & fortunato 2012, sci. rep.)
for j = 1:size(cipks,3)
    cicon(:,j) = fcn_consensus_communities(cipks(:,:,j),nreps,true); %  these represent the final clusters -- you'll get multiple resolutions (different number/size of clusters) that will be more/less-suited for behavioral analysis
end


%% save
% network_an = rmfield(network,'R');
netanal.data = data;
netanal.rho = R;
netanal.lag = lag;
netanal.ci = ci;
netanal.gamma = gammavals;
if size(cicon,2)>length(t)
    t(length(t)+1:size(cicon,2)) = t(end);
end
netanal.thresholds = t;
netanal.conCom = cicon;
s = zeros(1,size(cicon,2));
for comi = 1:size(cicon,2)
    s(comi) = length(unique(cicon(:,comi)));
end
netanal.comSizes = s;

% save(filename,'netanal','-append');