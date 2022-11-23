
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluating the Mutual Information (MI) between Schaefer ROIs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ParRunFiles store the presentation of each block for each par for each
% run
ParRunFiles=addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Data\Study2\SwitchStay\New\ParRunFiles');
% ROIFiles = the mean timeseries of activity extracted from filtered func
% data for each par, per run, per ROI
ROIFiles=addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Data\Study2\SwitchStay\New\SchaeferROIs');
% Add MI package
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\mi');
% Add ICC function
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\ICC');
% Add InvRankTrans
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Functions');
% For plots
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\ColorBrewer');
% For labels
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Data');
% For. . . a function?
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA')
%For Motion Outliers
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Subj22MotionOutliers')

num_blocks=8;
num_pars=27;
num_rois=100;
%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a standard timeser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The task structure is 16 s of rest, followed by
task (94 sec), then rest, and repeat. TR=2 sec,
94/2=47, so 47 vols per task block. and 16/2=8, so 8 vols per rest.
So alternating 8, 47, 8, 47, etc
I'll make rest = 0, and all other task blocks = 1:8
%}

Rest=ones(8,1)*0;
B1=ones(47,1)*1;
B2=ones(47,1)*2;
B3=ones(47,1)*3;
B4=ones(47,1)*4;
B5=ones(47,1)*5;
B6=ones(47,1)*6;
B7=ones(47,1)*7;
B8=ones(47,1)*8;

time=[Rest;B1;Rest;B2;Rest;B3;Rest;B4;Rest;B5;Rest;B6;Rest;B7;Rest;B8;Rest];
time(end:451,1)=0;

% Create a structure that has the normalized timeseries of activity per subj, per run, per ROI. The final table also has a column of what block each volume was in. 
SchaeferMaster={};

% Loop through subjs
for subj=1:num_pars
    if subj==11 || subj==17 || subj==18
        continue
    end

% Look through runs
for rr=1:3
if subj==9 && rr==3
    continue
end

% Load the .mat files specifying the order of blocks
cd ('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Data\Study2\SwitchStay\New\ParRunFiles')
sessions = dir(strcat(num2str(subj),'_Run*')); 
load(sessions(rr).name);

% In struct task, the blocks are not labelled with/ standard Blocks. 
% The next few lines will figure out which order this par's blocks were in
% for this run.
StructTask=struct2table(StructTask);

for ii=1:height(StructTask) 
    if StructTask.temporal(ii) == 1 && StructTask.perceptual(ii) == 1 && StructTask.spatial(ii) == 1
        parTimeSer(ii,1)=1;
    elseif StructTask.temporal(ii) == 1 && StructTask.perceptual(ii) == 1 && StructTask.spatial(ii) == 2
        parTimeSer(ii,1)=2;
    elseif StructTask.temporal(ii) == 1 && StructTask.perceptual(ii) == 2 && StructTask.spatial(ii) == 1
        parTimeSer(ii,1)=3;
    elseif StructTask.temporal(ii) == 1 && StructTask.perceptual(ii) == 2 && StructTask.spatial(ii) == 2
        parTimeSer(ii,1)=4;    
    elseif StructTask.temporal(ii) == 2 && StructTask.perceptual(ii) == 1 && StructTask.spatial(ii) == 1
        parTimeSer(ii,1)=5;    
    elseif StructTask.temporal(ii) == 2 && StructTask.perceptual(ii) == 1 && StructTask.spatial(ii) == 2
        parTimeSer(ii,1)=6;
    elseif StructTask.temporal(ii) == 2 && StructTask.perceptual(ii) == 2 && StructTask.spatial(ii) == 1
        parTimeSer(ii,1)=7;        
    else
        parTimeSer(ii,1)=8;
    end
end 

% Awesome! ParTimeSer now gives the order of blocks this par had on this
% run, and the var order makes it even simpler. 
% Let's load in the data- this gives me a TR x ROI table
order=unique(parTimeSer,'stable');

% Load the extracted timeseries - this gives me a TR x ROI table
cd ('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\Data\Study2\SwitchStay\New\SchaeferROIs\')
name=strcat('S',num2str(subj),'_R',num2str(rr),'.csv');
df=readtable(name);

% Create an array where each TR indicates which block it is using the order
% var we created. This is where rest blocks are added
partime=Rest;
for ii=1:8
    partime(height(partime)+1:height(partime)+47,1)=ones(47,1)*order(ii,1);
    partime(height(partime)+1:height(partime)+8)=Rest;
end

% Account for the HRF- shift the timeseries up 10 sec/5 frames, then add to the df
% This means we need 446 vols of the timeser
num_frames=6;
HRF=ones(num_frames,1)*0;
partimeSer=[HRF;partime(1:451-num_frames,1)];

% Zscore everything to have centre mean 0 and 1 SD
for ii=1:num_rois
    tmp=table2array(df(:,ii));
%     tmp=InvRankTrans(tmp);
    tmp=zscore(tmp);
    df(:,ii)=array2table(tmp);
end

% Add partimeSer to the master list - now we have a structure that has the
% timeseries AND block info per subj
df(:,width(df)+1)=array2table(partimeSer);
SchaeferMaster{subj,rr}=df;

clear parTimeSer StructTask df partime 
end
end

cd ('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\ROS_TaskSwitching\LEiDA\')

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute pairwise MI for all ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We're missing data for subjs 11, 17, and 18, so we skip them. 22 was
% shown to have really low reliability and a fair amount of motion, so we're removing them, too
%}

for subj=1:num_pars
if subj==11 || subj==17 || subj==18 || subj==22
    continue
end   

% I need to add a special calculation for par 9, or just exclude them, as I
% do here
for rr=1:3
if subj==9 && rr==3
    continue
end

tmp=table2array(SchaeferMaster{subj,rr});

for ii=1:num_rois
    for jj=1:num_rois
        MI{subj,rr}(ii,jj) = mutualinfo(tmp(:,ii),tmp(:,jj));
    end
end

end
end

SchaeferLabels=table2array(readtable('Schaefer100ROI17NetworkLabel.txt'));

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assess within subject reliability of each region's MI across the 3 runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Helpful link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4913118/ 
Helpful link: https://d1wqtxts1xzle7.cloudfront.net/25350178/mcgrawk1996a-with-cover-page-v2.pdf?Expires=1660642115&Signature=fzrXIpkclBHtFjBhws9swJ2eP1vGyvQYmmyTs4gMHTSeri5iDMTfndk3ArFwDd~aq5~qkDzTNz~ooiU6uSFpN7k0xn3Vu0s0UO7igibBS~Setz52zcm3DVA5vsql1mm0EJ2I1jo11do0J80vtsB4d--HmQObGnSpmkNppHPowthRDBv1~yCbhJF2kTLy5NStzbYaryPyNnc37vnEDib4HzMfJCTu4Fnq7KIP-BZz-qmbHhMbzCTT6RHj2QkUJzelNUsBBOzUS3rqY6NQXOowrKw-JGCACknx84pctAXCTrl-o92k4TCCvVUtaiHfcDJ~drDiKKsfS7zlfHiN6dJadQ__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA

There are three things you must answer before running an ICC:
Model: One-Way Random Effects, Two-Way Random Effects, or Two-Way Mixed
Effects. We have 2-Way random effects, because of the random effects of
subject and block order. We do not expect a mixed effect (interaction?),
but even so, the stats are the same, it's just interpretative
differences.
Type of Relationship: Consistency or Absolute Agreement. We want absolute
agreement, since we are not concerned about additive relationships, just
the prescence of differences. 
Unit: Single rater or the mean of rater. We will use mean rater,
because I think we'll average both within and across subjects in future
analyses. 

%}
for subj=1:num_pars 
if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
continue
end  

for nn=1:num_rois

r1=MI{subj,1}(:,nn);
r2=MI{subj,2}(:,nn);
r3=MI{subj,3}(:,nn);
rr=[r1,r2,r3];

% So to summarize - we want a 2-way, absolute, mean rater method, Case 2A in McGraw, or 'A-1' in the ICC function 
[WhnSubjRelR(subj,nn), ~, ~, ~, ~, ~, WhnSubjRelP(subj,nn)] = ICC(rr,'A-1');
%[S22(1,nn), ~, ~, ~, ~, ~, ~]=ICC(rr,'A-1');
end
end

% Delete the subjs we don't need - messes up the visualization
WhnSubjRelR([9,11,17,18,22],:)=[];
WhnSubjRelP([9,11,17,18,22],:)=[];

% Visualize the ICC in MI for all regions for all subjs
figure()
colormap(brewermap([],"YlGnBu"));
h=imagesc(WhnSubjRelR);
colorbar()
xticks(0:2:100);
xticklabels(SchaeferLabels);

forMean=reshape(WhnSubjRelR,[2200,1]);
mm=mean(forMean);
sstd=std(forMean);
disp(strcat('Mean WhnSubjRel=',num2str(mm),', std=',num2str(sstd)))

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assess between subject reliability of each region's MI across subjs for each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We want to see whetehr the MI for each Region across participants is
similar - so for each run we need to make a table where columns are the
MI for a subject for one region.

%}

for run=1:3
for nn=1:num_rois
for subj=1:num_pars

% Fill the columns with ones for the subjs we don't have data for
if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
   rr(:,subj)=ones(height(rr),1);
else   
    rr(:,subj)=MI{subj,run}(:,nn);
end

end
% Delete the cols for subjs we don't have data for 
rr(:,[9,11,17,18,22])=[];

% We want a 2-way, absolute, mean rater method, Case 2A in McGraw, or 'A-1' in the ICC function 
[BetSubjRelR(run,nn), ~, ~, ~, ~, ~, BetSubjRelP(run,nn)] = ICC(rr,'A-1');
end
end

% Visualize the ICC in MI for all regions for all subjs
figure()
colormap(brewermap([],"YlGnBu"));
h=imagesc(BetSubjRelR);
colorbar()
xticks(0:2:100);
yticks(0:2:100);
xticklabels(SchaeferLabels);

forMean=reshape(BetSubjRelR,[300,1]);
mm=mean(forMean);
sstd=std(forMean);
disp(strcat('Mean BetwnSubjRel=',num2str(mm),', std=',num2str(sstd)))

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up to assess the effect of switch dims on MI
% Now we're averaging across runs for each ROI, for each par
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Helpful link:  https://uk.mathworks.com/matlabcentral/answers/1629745-repeated-measure-anova-in-matlab 
Originally, we computed the miFC for each run, per participant. Now we're doing that, but for each block. So instead of one 100x100 miFC table per person, per run, we would have 8.
So, steps. (1) For each subj, I need to group the ROI timeseries by block across runs, so I can average them. (2) For each subj, for each ROI, average the ROI timeseries across runs. This should result in a TR x ROI table for each block. (3) Compute the pairwise MI to get the 100x100 miFC table per person, one table per block.
Then we can run a repeated measures ANOVA assessing whether there's an effect of any switch dims on the MI between ROIs
%}

CrossRuns={};

for block=1:8
for par=1:num_pars
    
if par==9 || par==11 || par==17 || par==18 || par==22    
    continue
else
    
% SchaeferMaster{subj,run)
% Identify where Block 1 is for each ROI
tmp1 = table2array(SchaeferMaster{par,1}(:,101))==block;
tmp2 = table2array(SchaeferMaster{par,2}(:,101))==block;
tmp3 = table2array(SchaeferMaster{par,3}(:,101))==block;

% Extract timeseries of activity during Block 1 for each ROI
b_r1=table2array(SchaeferMaster{par,1}(tmp1,1:100));
b_r2=table2array(SchaeferMaster{par,2}(tmp2,1:100));
b_r3=table2array(SchaeferMaster{par,3}(tmp3,1:100));

% Average each ROIs activity across runs
for ii=1:num_rois
    tmp=[b_r1(:,ii),b_r2(:,ii),b_r3(:,ii)];
    B(:,ii)=mean(tmp,2);
end

% Put it in a structure CrossRuns{subj,block}
CrossRuns{par,block}=B;
end
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the MI among ROIs for each participant, for each block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BlockMI={};
for block=1:8
for subj=1:num_pars
if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
    continue
else

tmp=CrossRuns{subj,block};

for ii=1:num_rois
    for jj=1:num_rois
        BlockMI{subj,block}(ii,jj) = mutualinfo(tmp(:,ii),tmp(:,jj));
    end
end
end
end
end

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the repeated measure anova to assess the effect of switch dims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To see the effect of Switch Dims, I could run one ANOVA per unique ROxROI 
to see if there's an effect of the switch dims on their MI. 
This would require having a subj x block table per ROI-ROI pair. 
%}

clear tmp
warning('off','all')
% [AT,SeqPredSig,SpaPredSig,PerPredSig] = RunAnovas(num_rois,num_pars,BlockMI);
% warning('on','all')

for rr1=1:num_rois
for rr2=1:num_rois
for subj=1:num_pars
if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
    continue
else

for block=1:8
ROIMI(subj,block)=array2table(BlockMI{subj,block}(rr1,rr2));
end
end
end

% Delete the pars you didn't include
ROIMI([9,11,17,18,22],:)=[];
ROIMI.Properties.VariableNames = {'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8'};

% Now create the effects design
withinDesign = table([1 1 1 1 2 2 2 2]',[1 2 1 2 1 2 1 2]',[1 1 2 2 1 1 2 2]','VariableNames',{'SeqPred','SpaPred','PerPred'});
withinDesign.SeqPred = categorical(withinDesign.SeqPred);
withinDesign.SpaPred = categorical(withinDesign.SpaPred);
withinDesign.PerPred = categorical(withinDesign.PerPred);

%Create the repeated measures ANOVA model adn display outputs
rm = fitrm(ROIMI, 'b1-b8 ~ 1', 'WithinDesign', withinDesign);
AT = ranova(rm, 'WithinModel', 'SeqPred*PerPred*SpaPred');
% disp(anovaTable(AT, 'Measurement'))

% Now I need to save the significant outputs - I'll just save the main
% effects per dim

% Computing effect size, partial Eta2, which is SS/(SS+SSErrow), where SS=Sum of
% Squares
% For seq pred
if table2array(AT(3,8))<0.05
    SeqPredSig(rr1,rr2)=table2array(AT(3,8));
    Eta=table2array(AT(3,1))/(table2array(AT(3,1))+table2array(AT(4,1)));
    SeqPredEta(rr1,rr2)=Eta;

% for spa pred
elseif table2array(AT(5,8))<0.05
    SpaPredSig(rr1,rr2)=table2array(AT(5,8));
    Eta=table2array(AT(5,1))/(table2array(AT(5,1))+table2array(AT(6,1)));
    SpaPredEta(rr1,rr2)=Eta;
    
% for per pred
elseif table2array(AT(7,8))<0.05
    PerPredSig(rr1,rr2)=table2array(AT(7,8));
    Eta=table2array(AT(7,1))/(table2array(AT(7,1))+table2array(AT(8,1)));
    PerPredEta(rr1,rr2)=Eta;
    
else   
end
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct for multiple comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SpaPredSig2=[];
SeqPredSig2=[];
PerPredSig2=[];
UseSpaPred=exist('SpaPredSig','var');

if UseSpaPred == 1
for ii=1:width(SpaPredSig)
    [~,~,SpaPredSig2(:,ii)] = fdr(SpaPredSig(:,ii));
end
end

for ii=1:width(SeqPredSig)
    [~,~,SeqPredSig2(:,ii)] = fdr(SeqPredSig(:,ii));
end

for ii=1:width(PerPredSig)
    [~,~,PerPredSig2(:,ii)] = fdr(PerPredSig(:,ii));
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PER PRED Run post hocs for main effects of dimension on sig miFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear tmp tstat Easy Hard
p=[];
PerPredMatrixx=[];
PerPredPostHoc=[];

EasyGreater=0;
EasyLower=0;

for ii=1:width(PerPredSig2)
for jj=1:ii
    if PerPredSig2(ii,jj) < 0.05 &&  PerPredSig2(ii,jj) > 0 % If this roi-roi connection is significant,
        
        for subj=1:num_pars
        if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
            continue
        else
        for block=1:8
        val=BlockMI{subj,block}(ii,jj);
        tmp(subj,block)=val;

        end
        end
        end

        % Delete the pars you didn't include
        tmp([9,11,17,18,22],:)=[];

        Easy=tmp(:,[1,2,5,6]);
        Hard=tmp(:,[3,4,7,8]);
        
        Easy=reshape(Easy,[height(Easy)*4,1]);
        Hard=reshape(Hard,[height(Hard)*4,1]);
        h=height(p)+1;
        [sigVal,~,stats]=signrank(Easy,Hard);
        zstat(h,1)=stats.zval;
        p(h,1)=sigVal;
        
        if mean(Easy) > mean(Hard)
            EasyGreater=EasyGreater+1;
        else
            EasyLower=EasyLower+1;
        end
        
        if sigVal < 0.05
            PerPredMatrixx(ii,jj)=sigVal;
            PerPredPostHoc(h,1)=sigVal;
            PerPredPostHoc(h,2)=zstat(h,1);
        else
            PerPredPostHoc(h,1)=sigVal;
            PerPredPostHoc(h,2)=zstat(h,1);       
        end
       
    end
    clear tmp
end
end
% PerPredPostHoc(:,1)=p(:,1);
% PerPredPostHoc(:,2)=zstat(:,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEQ PRED Run post hocs for main effects of dimension on sig miFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear tmp tstat Easy Hard
p=[];
SeqPredMatrixx=[];
SeqPredPostHoc=[];
totalSigcount=0;
EasyGreater=0;
EasyLower=0;

for ii=1:width(SeqPredSig2)
for jj=1:ii
    if SeqPredSig2(ii,jj) < 0.05 &&  SeqPredSig2(ii,jj) > 0 % If this roi-roi connection is significant,
        
        for subj=1:num_pars
        if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
            continue
        else
            
        for block=1:8
        val=BlockMI{subj,block}(ii,jj);
        tmp(subj,block)=val;

        end
        end
        end

        % Delete the pars you didn't include
        tmp([9,11,17,18,22],:)=[];
        
        Easy=tmp(:,[1,2,3,4]);
        Hard=tmp(:,[5,6,7,8]);
        Easy=reshape(Easy,[height(Easy)*4,1]);
        Hard=reshape(Hard,[height(Hard)*4,1]);
        h=height(p)+1;
        [p(h,1),~,stats]=signrank(Easy,Hard);
        zstat(h,1)=stats.zval;
        
        if mean(Easy) > mean(Hard)
            EasyGreater=EasyGreater+1;
        else
            EasyLower=EasyLower+1;
            disp(strcat('EasyLower for ROI',num2str(ii),'and',num2str(jj)))
        end
        
        if p(h,1) < 0.05
            totalSigcount=totalSigcount+1;
            SeqPredMatrixx(ii,jj)=p(h,1);
            SeqPredPostHoc(h,1)=p(h,1);
            SeqPredPostHoc(h,2)=zstat(h,1);
        else
            SeqPredPostHoc(h,1)=p(h,1);
            SeqPredPostHoc(h,2)=zstat(h,1);        
        end
%         forFig=[Easy,Hard]
%         figure()
%         hold on
%         boxplot(forFig)
        
    end
    clear tmp
end
end
% SeqPredPostHoc(:,1)=p(:,1);
% SeqPredPostHoc(:,2)=zstat(:,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPA PRED Run post hocs for main effects of dimension on sig miFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear tmp tstat Easy Hard
p=[];
SpaPredMatrixx=[];
if UseSpaPred == 1
for ii=1:width(SpaPredSig2)
for jj=1:ii
    if SpaPredSig2(ii,jj) < 0.05 &&  SpaPredSig2(ii,jj) > 0 % If this roi-roi connection is significant,
        
        for subj=1:num_pars
        if subj==9 || subj==11 || subj==17 || subj==18 || subj==22
            continue
        else
        for block=1:8
        val=BlockMI{subj,block}(ii,jj);
        tmp(subj,block)=val;

        end
        end
        end

        % Delete the pars you didn't include
        tmp([9,11,17,18,22],:)=[];
        
        Easy=tmp(:,[1,3,5,7]);
        Hard=tmp(:,[2,4,6,8]);
        Easy=reshape(Easy,[height(Easy)*4,1]);
        Hard=reshape(Hard,[height(Hard)*4,1]);
        h=height(p)+1;
        [p(h,1),~,stats]=signrank(Easy,Hard);
        zstat(h,1)=stats.zval;
        if p(h,1) < 0.05
            SpaPredMatrixx(ii,jj)=p(h,1);
            SpaPredPostHoc(h,1)=p(h,1);
            SpaPredPostHoc(h,2)=zstat(h,1);
        else
            SpaPredPostHoc(h,1)=nan;
            SpaPredPostHoc(h,2)=nan;
        end
        
    end
    clear tmp
end
end
% SpaPredPostHoc(:,1)=p(:,1);
% SpaPredPostHoc(:,2)=zstat(:,1);
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save effect sizes that survive post hocs for figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PerPredEta2=[];
% SeqPredEta2=[];
% SpaPredEta2=[];
% 
% for ii=1:height(PerPredMatrixx)
% for jj=1:width(PerPredMatrixx)
%     if PerPredMatrixx(ii,jj) < 0.05 && PerPredMatrixx(ii,jj) > 0
%         PerPredEta2(ii,jj) = PerPredEta(ii,jj);
%     else
%         PerPredEta2(ii,jj) = 0 ;
%     end
% end
% end
% 
% if UseSpaPred == 1
% for ii=1:height(SpaPredEta)
% for jj=1:width(SpaPredEta)
%     if SpaPredMatrixx(ii,jj) < 0.05
%         SpaPredEta2(ii,jj) = SpaPredEta(ii,jj);
%     else
%         SpaPredEta2(ii,jj) = 0 ;
%     end
% end
% end
% end
% 
% for ii=1:height(SeqPredMatrixx)
% for jj=1:width(SeqPredMatrixx)
%     if SeqPredMatrixx(ii,jj) < 0.05 && SeqPredMatrixx(ii,jj) > 0
%         SeqPredEta2(ii,jj) = SeqPredEta(ii,jj);
%     else
%         SeqPredEta2(ii,jj) = 0 ;
%     end
% end
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save effect sizes that survive fdr correction for figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PerPredEta2=[];
SeqPredEta2=[];
SpaPredEta2=[];

for ii=1:height(PerPredSig2)
for jj=1:width(PerPredSig2)
    if PerPredSig2(ii,jj) < 0.05 && PerPredSig2(ii,jj) > 0
        PerPredEta2(ii,jj) = PerPredEta(ii,jj);
    else
        PerPredEta2(ii,jj) = 0 ;
    end
end
end

for ii=1:height(SeqPredSig2)
for jj=1:width(SeqPredSig2)
    if SeqPredSig2(ii,jj) < 0.05 && SeqPredSig2(ii,jj) > 0
        SeqPredEta2(ii,jj) = SeqPredEta(ii,jj);
    else
        SeqPredEta2(ii,jj) = 0 ;
    end
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The arrays aren't a full 100 ROIs wide- I'll fix that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if UseSpaPred == 1
SpaPredEta3=SpaPredEta2;
wneeded=num_rois-width(SpaPredEta2);
hneeded=num_rois-height(SpaPredEta2);
SpaPredEta3(:,width(SpaPredEta3)+1:num_rois)=ones(height(SpaPredEta3),wneeded)*0';
SpaPredEta3(height(SpaPredEta3)+1:num_rois,:)=ones(hneeded,num_rois)*0;
end

SeqPredEta3=SeqPredEta2;
wneeded=num_rois-width(SeqPredEta2);
hneeded=num_rois-height(SeqPredEta2);
SeqPredEta3(:,width(SeqPredEta3)+1:num_rois)=ones(height(SeqPredEta3),wneeded)*0';
SeqPredEta3(height(SeqPredEta3)+1:num_rois,:)=ones(hneeded,num_rois)*0;

PerPredEta3=PerPredEta2;
wneeded=num_rois-width(PerPredEta2);
hneeded=num_rois-height(PerPredEta2);
PerPredEta3(:,width(PerPredEta3)+1:num_rois)=ones(height(PerPredEta3),wneeded)*0';
PerPredEta3(height(PerPredEta3)+1:num_rois,:)=ones(hneeded,num_rois)*0;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To make circgraph with the proper alignment, I've got to reorganize the
% data and the labels. 
SchaeferLabels2(1:50,1)=SchaeferLabels(1:50,1);
Flipped=flip(SchaeferLabels(51:100,1));
SchaeferLabels2(51:100,1)=Flipped;
Flipped=[];

% Cool, now just gotta flip results the same way
SeqPredEta4=SeqPredEta3;
SeqPredEta4(51:end,1:50)=flip(SeqPredEta4(51:end,1:50));
SeqPredEta4(1:50,51:end)=flip(SeqPredEta4(1:50,51:end),2);
SeqPredEta4(51:end,51:end) = flip(flip(SeqPredEta4(51:end,51:end),2));

PerPredEta4=PerPredEta3;
PerPredEta4(51:end,1:50)=flip(PerPredEta4(51:end,1:50));
PerPredEta4(1:50,51:end)=flip(PerPredEta4(1:50,51:end),2);
PerPredEta4(51:end,51:end) = flip(flip(PerPredEta4(51:end,51:end),2));


% Now I have to make a color map where the first 50 vals map
C=linspecer(17);
load Clist
for ii=1:height(Clist)
   colr(1,:)= C(Clist(ii,1),:);  
   Clist2(ii,:)=colr ;
end

% Plot some sheeeet
figure()
circularGraph(SeqPredEta4,'Label',SchaeferLabels2,'Colormap',Clist2)

if UseSpaPred == 1
figure()
circularGraph(SpaPredEta4,'Label',SchaeferLabels)
end

figure()
circularGraph(PerPredEta4,'Label',SchaeferLabels2,'Colormap',Clist2)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create tables for p and effect size per miFC connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SigTableSeqPredROIs={};
SigTablePerPredROIs={};
SigTableSpaPredROIs={};
SigTableSeqPred=[];
SigTablePerPred=[];
SigTableSpaPred=[];

strang={};

for ii=1:height(SeqPredSig2)
for jj=1:ii
   if SeqPredSig2(ii,jj) < 0.05 && SeqPredSig2(ii,jj) > 0
       h=height(SigTableSeqPred)+1;
       SigTableSeqPredROIs(h,1)=SchaeferLabels(ii);
       SigTableSeqPredROIs(h,2)=SchaeferLabels(jj);

       SigTableSeqPred(h,1)=SeqPredSig2(ii,jj);
       SigTableSeqPred(h,2)=SeqPredEta(ii,jj);
       SigTableSeqPred(h,3)=ii;
       SigTableSeqPred(h,4)=jj;   
   end
end
end
   
   
if UseSpaPred == 1
for ii=1:height(SpaPredSig2)
for jj=1:ii
   if SpaPredSig2(ii,jj) < 0.05 && SpaPredSig2(ii,jj) > 0
       h=height(SigTableSpaPred)+1;
       SigTableSpaPredROIs(h,1)=SchaeferLabels(ii);
       SigTableSpaPredROIs(h,2)=SchaeferLabels(jj);
       SigTableSpaPred(h,1)=SpaPredSig2(ii,jj);
       SigTableSpaPred(h,2)=SpaPredEta(ii,jj);
   end   
end
end
end

for ii=1:height(PerPredSig2)
for jj=1:ii
   if PerPredSig2(ii,jj) < 0.05 && PerPredSig2(ii,jj) > 0
       h=height(SigTablePerPred)+1;
       SigTablePerPredROIs(h,1)=SchaeferLabels(ii);
       SigTablePerPredROIs(h,2)=SchaeferLabels(jj);
       SigTablePerPred(h,1)=PerPredSig2(ii,jj);
       SigTablePerPred(h,2)=PerPredEta(ii,jj);
       SigTablePerPred(h,3)=ii;
       SigTablePerPred(h,4)=jj;
       %        disp(strcat('miFC',num2str(ii),',',num2str(jj),',',num2str(PerPredSig2(ii,jj)),',',num2str(PerPredEta3(ii,jj))));

   end
end
end

save('EverythingFromMI_HRF6_v2.mat','-v7.3')

