function GetBehavior(taskName,HeaderForCSV,TaskFiles)

%{

This script extracts the median reaction time and mean accuracy per participants who's
imaging data is being used. 

%}



% Get Subject names 
names={};
for k = 1:length(TaskFiles)
    names{k,1}=TaskFiles(k).name;
    names{k,1} = erase(names{k,1},".csv");   
    names{k,1} = erase(names{k,1},HeaderForCSV); 
 
end 

% Load table with behavioral data
TT=readtable("HCPBehavior.csv");

taskName=char(taskName);
MedRT=[];
Acc=[];
index=0;
for ii=1:height(TT)
    for jj=1:height(names)
        if TT.Subject(ii) == str2double(names{jj})
            index=index+1;
            
            switch taskName
                case "WM"
                    MedRT(index,1)=TT.WM_Task_Median_RT(ii);   % be sure to edit the name of the task as needed. 
                    Acc(index,1)=TT.WM_Task_Acc(ii);           % be sure to edit the name of the task as needed. 

                case "Lang"
                    MedRT(index,1)=TT.Language_Task_Median_RT(ii);   % be sure to edit the name of the task as needed. 
                    Acc(index,1)=TT.Language_Task_Acc(ii);           % be sure to edit the name of the task as needed. 

                case "Relat"
                    MedRT(index,1)=TT.Relational_Task_Median_RT(ii);   % be sure to edit the name of the task as needed. 
                    Acc(index,1)=TT.Relational_Task_Acc(ii);           % be sure to edit the name of the task as needed. 

                case "Gamb"
                    MedRT(index,1)=TT.Gambling_Task_Perc_Larger(ii);   % be sure to edit the name of the task as needed. 
                    Acc(index,1)=TT.Relational_Task_Acc(ii);           % be sure to edit the name of the task as needed. 
                    MedRT_small(index,1)=TT.Gambling_Task_Median_RT_Smaller(ii);   % be sure to edit the name of the task as needed. 
                    Acc_small(index,1)=TT.Gambling_Task_Perc_Smaller(ii);           % be sure to edit the name of the task as needed. 

                otherwise
                    disp('Please edit the function to include this task name!')

            end    
            
        else
         continue
        end
    end
end

loc=pwd;
NAME=strcat(loc,'\Outputs\','TaskPerformance_',taskName,'.mat');

switch taskName
    case "WM"
        save(NAME, 'MedRT', 'Acc')

    case "Lang"
        save(NAME, 'MedRT', 'Acc')

    case "Relat"
        save(NAME, 'MedRT', 'Acc')

    case "Gamb"
        save(NAME, 'MedRT', 'Acc','MedRT_small','Acc_small')
     
    otherwise
        disp('Please edit the function to include this task name!')

end

end

