function MM = markovmodel(X, N)
%% MARKOVMODEL 
% Return an Nth order Markov model representation of the state sequence
% X where X is a row vector of length L or a matrix (L row * M columns) 
% by Greg
if nargin < 2
    N = 1;
end

[C1,~,S1] = unique(X,'rows');
M = [];
for(t=1+N:length(S1))
    M = [M; S1(t-N:t)'];
end

[C2,~,S2] = unique(M(:,1:end-1),'rows');
T = sortrows([S2 M(:,end)]);

MM = [];
for i=1:T(end,1)
    frtab = tabulate(T(T(:,1) == i,2)); frtab(frtab(:,2) == 0,:) = [];
    MM = vertcat(MM, ...
        [ reshape(C1(C2(repmat(i,size(frtab,1),1),:)),size(frtab,1),[]), ...
        (C1(frtab(:,1))),...
        frtab(:,3)/100 ]); 
end
end