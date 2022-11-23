function GetReadyForLEiDA(taskName,TaskRest)

% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp('%%                     Preparing data for MI                       %%')
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

n = length(TaskRest); %this is the number of runs we're looking at

%%% REST
R = [];
for k = 1:n
    R(:,:,k) = load(TaskRest(k).name);  
end 

% Because my tables are flipped (rows=frame, column=ROI), this ensures the 
% tables so now it's in the proper format: each row is a ROI, each column 
% is a frame in the timeseries. This also ensures the rest timeseries is 
% the same length as task

cutoff= height(squeeze(R(:,:,1)));

df={};

for k=1:n
    
    % Flip rest
    xx = array2table(R(:,:,k));
    xx = rows2vars(xx) ;
    xx(:,1) = [];
    df{k,1}=table2array(xx);
    
end

% save DataForLEiDA_RestAndLang.mat
loc=pwd;
NAME=strcat(loc,'\Outputs\','DataForMI_',taskName,'.mat');
save(NAME, 'df', 'R')

end

