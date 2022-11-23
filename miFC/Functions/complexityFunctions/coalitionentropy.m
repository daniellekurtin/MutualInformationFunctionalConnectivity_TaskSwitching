function [ Hc, coals] = coalitionentropy(ths, Cth)
if nargin < 2
    Cth = 0.8;
end
% For data Xt, consisting of channels Xi,t, i = 1, ?, n, we consider two channels to be in synchrony at time t if the absolute value of the difference between their instantaneous Hilbert phases is less than 0.8 radians (approximately 45 degrees). Then we define coalition time-series  by  taking the value 1 if channels i and j are synchronised at time t and taking the value 0 otherwise   

% Coalition entropy measures the variety of metastable states entered by
% a system of oscillators.  As with the synchrony
% metric, we calculated the phase of each oscillator at each time point
% t using a Hilbert transform. 
%We then performed clustering at each time
% point by picking the two most synchronous oscillators/coalitions using
% the first equation defined for the synchrony metric. Once a pair was %
%identified they were joined to form a new coalition and the new coalitions
% mean complex exponential phase was calculated for use in the future most
% synchronous pair selection process. A threshold of 0.05 from full synchrony
% was used to limit the cluster merging. The process was repeated until no
% oscillators/coalitions fell within the threshold to allow merging into
% a new coalition.
T = size(ths,1);
N = size(ths,2);

coals = zeros(T,N);
for t=1:T
    coals(t,:) = calccoalitions(ths(t,:), Cth);
end

[C,IA,IC] = unique(coals,'rows');
f = tabulate(IC);
p = f(:,3)/100;

S = T;
Hc = (-sum(p.*log2(p)))/log2(S);
end
