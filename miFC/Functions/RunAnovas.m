function [AT,SeqPredSig,SpaPredSig,PerPredSig] = RunAnovas(num_rois,num_pars,BlockMI)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
clear tmp
warning('off','all')

for rr1=1:num_rois
for rr2=1:num_rois
for subj=1:num_pars
if subj==9 || subj==11 || subj==17 || subj==18
    continue
else

for block=1:8
ROIMI(subj,block)=array2table(BlockMI{subj,block}(rr1,rr2));
end
end
end

% Delete the pars you didn't include
ROIMI([9,11,17,18],:)=[];
ROIMI.Properties.VariableNames = {'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8'};

% Now create the effects design
withinDesign = table([1 1 1 1 2 2 2 2]',[1 2 1 2 1 2 1 2]',[1 1 2 2 1 1 2 2]','VariableNames',{'SeqPred','SpaPred','PerPred'});
withinDesign.SeqPred = categorical(withinDesign.SeqPred);
withinDesign.SpaPred = categorical(withinDesign.SpaPred);
withinDesign.PerPred = categorical(withinDesign.PerPred);

%Create the repeated measures ANOVA model adn display outputs
rm = fitrm(ROIMI, 'b1-b8 ~ 1', 'WithinDesign', withinDesign);
AT = ranova(rm, 'WithinModel', 'SeqPred*PerPred*SpaPred');
end

if table2array(AT(3,8))<0.045
    SeqPredSig(rr1,rr2)=table2array(AT(3,8));

% for spa pred
elseif table2array(AT(5,8))<0.045
    SpaPredSig(rr1,rr2)=table2array(AT(5,8));

% for per pred
elseif table2array(AT(7,8))<0.045
    PerPredSig(rr1,rr2)=table2array(AT(7,8));
else   
end
warning('on','all')
end
end