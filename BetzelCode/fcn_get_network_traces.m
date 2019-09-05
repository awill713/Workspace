function [Rho,Lag] = fcn_get_network_traces(data)
% FCN_GET_NETWORK       generate lagged correlation matrix
%
%   [Rho,Lag] = FCN_GET_NETWORK(data) takes a data structure and generates
%         a correlation matrix, Rho.
%
%   Inputs:
%   
%   data,       structure with the following fields:
%               data.ts = [neuron x time] matrix of traces
%               data.proc.filt.filttype = filter type ('low', 'high',
%                   'bandpas', or 'stop')
%               data.proc.filt.freqrange = frequences for filter
%               data.proc.diff = true/false to difference time series
%               data.proc.gsr = true/false to perform global signal regions
%               data.net.maxlag_seconds = maximum lag (in seconds) over
%                   which we compute correlations.
%               data.fr = sampling frequency in Hertz
%
%   Outputs:
%
%   Rho,        neuron x neuron correlation matrix
%   Lag,        neuron x neuron matrix of time lags
%
%   Richard Betzel, UPenn 2017
%
%   Note: this procedure is similar to what is described in Muldoon et al
%   2013, PNAS. Diff: we do not convolve spike trains with kernels because
%   calcium imaging is already slow/continuous signal.
%    
%   A couple additional/possible processing steps
%% get dimensions of data
[N,NT] = size(data.ts);             % number of nodes/neurons and time points
%% global signal regression
if data.proc.gsr
    mu = nanmean(data.ts);          % orthogonalize neuronal time series wrt global signal
    for i = 1:N
        [~,~,data.ts(i,:)] = regress(data.ts(i,:)',[ones(NT,1),mu']);
    end
end
%% filtering
if any(data.proc.filt.filttype > 0)
    data.ts = ...                        % filter the traces
        buttfilt(...
        data.ts,...                     % orig time series
        data.proc.filt.freqrange,...    % frequencies for filter
        data.fr,...                     % sampling frequency
        data.proc.filt.filttype,...     % filter type: bandpass/low/high/stop
        data.proc.filt.filtorder);     % filter order
end
%% differencing
if data.proc.diff
    data.ts = diff(data.ts,[],2);          % difference filtered time series -- helps remove serial/auto correlation but makes interpretation more challenging
    NT = NT - 1;
end
%% zero-lag covariance
data.ts = bsxfun(@minus,data.ts,mean(data.ts,2));       % zero-mean time series
maxlag = ceil(data.fr*data.net.maxlag_seconds);         % maximum lab in samples
cov_std = zeros(N,1);                                   % allocate memory for vector of zero-lag stdv
for i = 1:N                                             % loop over nodes
    cov_std(i) = xcov(data.ts(i,:),data.ts(i,:),0);       % compute variance and take sqrt for standard deviation
end
cov_std = sqrt(cov_std);
normalization = cov_std*cov_std';                       % we'll use this to normalize the lagged covariance matrix
%% pairwise lagged covariance
Rho = zeros(N);                                         % allocate memory for correlation matrix
Lag = zeros(N);                                         % allocate memory for lag matrix
for i = 1:(N - 1)                                       % loop i nodes
    for j = (i + 1):N                                   % loop j nodes
        r = xcov(data.ts(i,:),data.ts(j,:),maxlag);       % lagged covariance
        [~,Lag(i,j)] = max(abs(r));                     % find lag at which max absolute cov. occurs
        Rho(i,j) = r(Lag(i,j));                         % store this value
    end
    imagesc(Rho./normalization); drawnow;
end
%% final processing
Rho = (Rho + Rho')/2;                                   % fill lower triangle
Rho = Rho./normalization;                               % cov -> corr
Rho = 0.5*log((1 + Rho)./(1 - Rho));                    % fisher r -> z transform
Rho(1:(N + 1):end) = 0;                                 % make diagonal zero