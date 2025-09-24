%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Plot k-range ---
% mCAPs script based on the work of Farnaz Delavari (mCAP-main)
% plot tested k range
%
%
% Created by : Camille Serquet 
% Creation : 05.2025
% Last modification : -
% MATLAB version : R2022b
% SPM version : 12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*** START ***')

%% SETTING ENVIRONMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderPath = fullfile('./in_outputs/findK');
outPath = fullfile('./in_outputs');
CPUpath = fullfile('./in_outputs/CPUinfos.xlsx');
%% EXTRACT *.MAT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foldersList = struct2table(dir(fullfile(folderPath, 'TT*'))).name; %get list of files in folder

sumTbl = zeros(size(foldersList,1),2);

for i = 1:size(foldersList,1)
    matFile = dir(fullfile(folderPath, foldersList{i}, '*.mat')); %get mat file in subfolder
    m = matfile(fullfile(folderPath, foldersList{i}, matFile.name)); %load data matrix
    
    sumTbl(i,2) = m.meanMaxSumD; %get mean value of the max sum of distances
    sumTbl(i,1) = str2double(extractBefore(extractAfter(matFile.name, 'K'), '.mat')); %extract k number
       
end

sumTbl = sortrows(sumTbl,1);
k = sumTbl(:,1);
meanMaxSumD = sumTbl(:,2);

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('WindowState','maximized');
s1 = axes(f);
hold(s1,'on')
% s2 = axes(f);
% hold(s2, 'on')
% s3 = axes(f);
% hold(s3, 'on')

%%----------k vs meanMaxSumD-----------------------------------------------
%subplot(1,3,1,s1)
plot(s1,k,meanMaxSumD, 'LineWidth', 2, 'Color', '#4E95D9','Marker','diamond','MarkerSize',15,'MarkerFaceColor','#4E95D9')

plot(s1, [6 6], [0 meanMaxSumD(4)], 'LineStyle', ':', 'LineWidth', 4, 'Color', '#E97132', 'Marker','diamond','MarkerSize',16,'MarkerFaceColor','#E97132')

title(s1, "Concatenated Means of Maximum Sum of Distances")

xlabel(s1, "Tested k (-)")
xticks(s1, min(k):1:max(k))
xlim(s1, [min(k) max(k)])

ylabel(s1, "Mean of Maximum Sum of Distances (-)")
yticks(s1, 1000:500:3000)
ylim(s1, [1000 3000])

grid(s1,'on')
set(s1,'FontSize',16)
set(s1, 'Color','#DCEAF7')

% %% [OPTIONAL]
% %%-----------k vs CPUefficiency--------------------------------------------
% 
% CPUinfo = readtable(CPUpath);
% k = CPUinfo.k;                                                              %k range
% mem  = CPUinfo.memory;                                              %memory utilized                                                  %time to compute
% t = hours(duration(CPUinfo.time,"InputFormat","dd:hh:mm:ss","Format","h"));            %convert time in hours
% subplot(1,3,2,s2)
% plot(s2,k,mem, 'LineWidth', 2, 'Color', '#4E95D9','Marker','diamond','MarkerSize',15,'MarkerFaceColor','#4E95D9')
% 
% plot(s2, [6 6], [0 mem(4)], 'LineStyle', ':', 'LineWidth', 4, 'Color', '#E97132', 'Marker','diamond','MarkerSize',16,'MarkerFaceColor','#E97132')
% 
% title(s2, "CPU Memory Used")
% 
% xlabel(s2, "Tested k (-)")
% xticks(s2, min(k):1:max(k))
% xlim(s2, [min(k) max(k)])
% 
% ylabel(s2, "CPU Memory (GB)")
% yticks(s2, 70:1:80)
% ylim(s2, [70 80])
% 
% grid(s2,'on')
% set(s2,'FontSize',16)
% set(s2, 'Color','#DCEAF7')
% 
% %%-----------k vs CPUefficiency--------------------------------------------
% subplot(1,3,3,s3)
% plot(s3,k,t, 'LineWidth', 2, 'Color', '#4E95D9','Marker','diamond','MarkerSize',15,'MarkerFaceColor','#4E95D9')
% 
% plot(s3, [6 6], [0 t(4)], 'LineStyle', ':', 'LineWidth', 4, 'Color', '#E97132', 'Marker','diamond','MarkerSize',16,'MarkerFaceColor','#E97132')
% 
% title(s3, "Time needed to compute K")
% 
% xlabel(s3, "Tested k (-)")
% xticks(s3, min(k):1:max(k))
% xlim(s3, [min(k) max(k)])
% 
% ylabel(s3, "Time (h)")
% yticks(s3, 5:5:20)
% ylim(s3, [5 20])
% 
% grid(s3,'on')
% set(s3,'FontSize',16)
% set(s3, 'Color','#DCEAF7')