function W = graph2net(sA,e,d)

% function W = graph2net(sA,e,d)
%
% sA = nxn binary adjacency matrix for a directed graph
%      -> can also have t values in [0,1] for interpolation b/w graphs
% e = epsilon value to control synaptic strengths (default: e = 0.25)
%     -> can be a vector of length n, for varied eps_i's, or a scalar
% d = delta value for inhibitory weights (default: d = 2*e)
%     -> can be a vector of length n, for varied delta_i's, or a scalar
%
% create matrix W for a threshold-linear network from directed graph
% -> specifically, this creates W for a CTLN or generalized CTLN
%
% modified june 22, 2015
% updated dec 30, 2020 to allow generalized sA matrices for interpolation
% updated apr 23, 2022 to fix bug: A was unassigned, now fixed
% updated aug 10, 2023 to remove A since unnecessary
% updated mar 29, 2025 to allow eps and delta (e,d) to be vectors,
% as in generalized CTLNs (gCTLNs)
% comments updated mar 30, 2025

n = size(sA,1);

if nargin < 2 || isempty(e)
    e = 0.25*ones(1,n);
end

if nargin < 3 || isempty(d)
    d = 0.5*ones(1,n);
end

% if e,d are not long enough, just use first entry of each
% to make them into constant vectors (this is for backwards compatibility)
if length(e)<n
    e = e(1)*ones(1,n);
end

if length(d)<n
    d = d(1)*ones(1,n);
end

% create matrix W from sA, allow e(j), d(j) to depend on source node j
% note that entries can be t in [0,1] to interpolate between graphs
W = zeros(size(sA));
for i=1:n
    for j=1:n
        t = sA(i,j);
        W(i,j) = t*(-1+e(j)) + (1-t)*(-1-d(j));
    end
    W(i,i) = 0; % make sure diagonal is 0
end