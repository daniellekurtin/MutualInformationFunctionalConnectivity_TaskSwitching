function [ Hroi, coalsroi] = sce_danielle(phases, thresh)
% SCE_DANIELLE Rough function to calculate the synchrony coalition entropy
% for each time series. NOTE THERE IS ANOTHER APPPROACH (SEE OTHER CODE)
% THAT CALCULATES A SINGLE 'GLOBAL' SCE
%
% INPUTS:
% phases = region x time ROI phases
% thresh = threshold that defines whether two phases are 'locked'
% i.e. For data Xt, consisting of channels Xi,t, i = 1, ?, n, we
% consider two channels to be in synchrony at time t if the absolute
% value of the difference between their instantaneous Hilbert phases
% is less than some threshold e.g 0.8
Nroi  = size(phases, 1);
Nvols = size(phases, 2); % the time dimension

% for each ROI, at each time point, work out whether the phase of the ROI
% is 'locked'/in phase with each other ROI

% you need to define the function that tells you if two phase are 'locked'
% I can't quite remember how to do this but you have to take into account
% that they wrap around the cirlce etc.
% areinsync = @(a,b) abs(a-b) < thresh; % THIS IS NOT CORRECT!!!!

areinsync = @(a,b) cos(a-b) > thresh; % This bounds it between 0 and 1

coalsroi = {};
for i=1:Nroi
    for t=1:Nvols
        % I'm going to do this very laboriously
        for k=1:Nroi
            coalsroi{i}(k,t) = areinsync(phases(i,t), phases(k,t));
        end
    end
end

% now work out the entropy values
Hroi = nan(Nroi,1);

for i=1:Nroi
    % this is just a way to turn our roi x time matrices
    % into a 1 x time list of 'states'
    [~,~,IC] = unique(coalsroi{i}','rows');
    % now we just tally up the frequency of each state
    % and turn it into a probability
    f = tabulate(IC);
    p = f(:,3)/100;
    % and entropy from that
    Hroi(i) = (-sum(p.*log2(p)));
end
end
