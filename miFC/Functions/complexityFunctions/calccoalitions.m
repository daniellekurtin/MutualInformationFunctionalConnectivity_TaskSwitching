function coals = calccoalitions(ths, thr)
N = length(ths);
coals = 1:N;

while(true)
    phis = zeros(N,N);
    
    % find the two most synchronous oscillators/coalitions
    for i=1:N
        thsi = angle(calcphi(ths(coals == i)'));
        for j=i+1:N
            phis(i,j) = ...
                calcphi([thsi, angle(calcphi(ths(coals == j)'))]');
        end
    end
    
    absphis   = abs(phis);
    maxabsphi = max(absphis(:));
    
    [i,j]     = find(absphis == maxabsphi, 1);
    
    if(i == j | maxabsphi < thr)
        break;
    end
    
    % once a pair identified, join them to form a new coalition
    coals(coals == j) = i; % all those in j take on  identity of i
    

end

end


% Coalition entropy measures the variety of metastable states entered by
% a system of oscillators.  As with the synchrony
% metric, we calculated the phase of each oscillator at each time point
% t using a Hilbert transform.
% We then performed clustering at each time
% point by picking the two most synchronous oscillators/coalitions using
% the first equation defined for the synchrony metric. Once a pair was %
% identified they were joined to form a new coalition and the new coalitions
% mean complex exponential phase was calculated for use in the future most
% synchronous pair selection process. A threshold of 0.05 from full synchrony
% was used to limit the cluster merging. The process was repeated until no
% oscillators/coalitions fell within the threshold to allow merging into
% a new coalition.
