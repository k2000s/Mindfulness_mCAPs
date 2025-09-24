%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Extract CAPs ---
% mCAPs script based on the work of Farnaz Delavari (mCAP-main)
% extracting CAPs nifti, seed nifti, and occurences of CAPs
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
spmPath = fullfile('./spm12');              %[EDIT]: path to SPM toolbox
dataPath = fullfile('./in_outputs/FINAL_mCAP_K6.mat');          %[EDIT]: path to data.mat
outPath = fullfile('./in_outputs');             %[EDIT]: path to folder to store data
tmpPath = fullfile('E:/databaseV2/func/wrfsub-s001_ses-01_task-rest_bold_.nii,1');
maskPath = fullfile('./in_outputs/r_brainMask.nii');
databasePath = fullfile('E:/databaseV2/func');
datamatPath = fullfile('./in_outputs/dataV2.mat');
slPath = '.\in_outputs\SubjectList.xlsx';
addpath(genpath(spmPath));                                                  %add path of SPM toolbox

if ~isfolder(outPath);  mkdir(outPath); end                                 %if output folder doesn't exist, create it

%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------load template for header----------------------------------------
tmpHeader = spm_vol(tmpPath); %load template header
tmpDim = tmpHeader.dim; %dimension of 3D volumes
%----------load brain mask-------------------------------------------------
brainMask = spm_read_vols(spm_vol(maskPath)); %read brain mask
brainMask(isnan(brainMask)) = 0;
brainMask = logical(brainMask); %convert to logical mask
%----------load seedmask---------------------------------------------------
m = matfile(datamatPath);
seedm = m.seedmask;
seedmask = zeros(size(brainMask));
seedmask(brainMask) = seedm;
%----------load output data matrix for selected k--------------------------
load(dataPath)

%% CREATE CAPS NIFTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it :  iterations made by kmeans until convergence
cap = M.CAPnorm{it-1};                                                      %load last CAP matrix
idxPos = find(brainMask);                                                   %get positive index in brain mask
for i = 1:size(cap,2)                                                       %for all CAPs    
    c = cap(:,i);                                                           %get i-th CAP
    capi = nan(size(brainMask));                                            %initiate i-th CAP volum
    for k = 1:numel(idxPos)
        capi(idxPos(k)) = c(k);
    end
        
    hdr = tmpHeader;                                                        %header of i-th CAP
    hdr.dt = [16 0];
    hdr.fname = fullfile(outPath,['mCAP_',char(num2str(i)),'.nii']);        %name of i-th CAP
       
    spm_write_vol(hdr, capi)                                                %write i-th CAP as nifti file
end

%% CREATE SEEDMAPS NIFTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parcel = 90;                                                                %parcel percentage
for i = 1:size(cap,2)                                                       %for all CAPs    
    c = cap(:,i);                                                           %get i-th CAP
    capi = zeros(size(brainMask));                                          %initiate i-th CAP volum
    for k = 1:numel(idxPos)                         
        capi(idxPos(k)) = c(k);
    end
    capi = capi.*seedmask;
    seedi = abs(capi(seedmask==1));
    thr = abs(capi)>=prctile(seedi,parcel);
    capSeed = capi.*thr;

    hdr = tmpHeader;                                                        %header of i-th CAP
    hdr.dt = [16 0];                                                        %range of grey scale
    hdr.fname = fullfile(outPath,['seedParcel90_',char(num2str(i)),'.nii']);%name of i-th CAP
       
    spm_write_vol(hdr, capSeed)                                             %write i-th CAP as nifti file
end



%% CALCULATE OCCURENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = M.indAll{it-1};                                                       %get the indices struct
kept = ind.kept.active;                                                     %get kept active voxel list
keptIdx = zeros(size(kept));                                                %initiate kept index with double values
idx = M.idx{it-1};                                                          %get cap number (to which caps it goes to)
keptIdx(kept==1) = idx;                                                     %put cap number on active indices

occRaw = zeros(size(cap,2), size(keptIdx,2));                               %initiate occurences table
for i = 1 : size(cap,2)                                                     %for each CAP
   occRaw(i,:) = sum(keptIdx==i);                                           %sum up the nb of occurences
end

%% CREATE CAPS XLSX FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------get subject & sessions number-----------------------------------
ls = dir(fullfile(databasePath,'w*'));                                      %get all func files
ls = struct2table(ls).name;                                                 %get all func files names

sub = nan(1,size(ls,1));                                                    %initiate subject nb vector
ses = nan(1,size(ls,1));                                                    %initiate session nb vector

for k = 1:size(ls,1)                                                        %for all file names
    %----------get subject number------------------------------------------
    subB = extractBefore(ls{k},'_');
    sub(k) = str2double(extractAfter(subB, '-s'));
    %----------get session number------------------------------------------
    sesA = extractAfter(ls{k}, 'ses-');
    ses(k) = str2double(extractBefore(sesA, '_'));
end

%-----------get group of each subjects-------------------------------------
slTable = readtable(slPath,"Sheet",'Session_timing','Range','A:B');
subjectList = table2array(unique(slTable));                                 %convert table to array and only keep unique CODE variables
subjectName = str2double(extractAfter(subjectList(:,1),'S'));
subjectCat = categorical(subjectList(:,2));                                
subCat = categorical(nan(size(ls, 1),1));

for i = 1:size(sub,2)
    for j = 1:size(subjectName, 1)
        if sub(i) == subjectName(j)
            subCat(i) = subjectCat(j);
        end
    end
end


CAP1_raw = occRaw(1,:)';     CAP2_raw = occRaw(2,:)';     CAP3_raw = occRaw(3,:)';
CAP4_raw = occRaw(4,:)';     CAP5_raw = occRaw(5,:)';     CAP6_raw = occRaw(6,:)';

%----------normalized CAPs for each ses------------------------------------
CAP1 = zeros(size(occRaw,2), 1);    CAP2 = zeros(size(occRaw,2), 1);
CAP3 = zeros(size(occRaw,2), 1);    CAP4 = zeros(size(occRaw,2), 1);
CAP5 = zeros(size(occRaw,2), 1);    CAP6 = zeros(size(occRaw,2), 1);
for k = 1:size(occRaw,2)
    subSum = sum(occRaw(1:6,k));
    CAP1(k) = CAP1_raw(k)*100/subSum;      CAP2(k) = CAP2_raw(k)*100/subSum;
    CAP3(k) = CAP3_raw(k)*100/subSum;      CAP4(k) = CAP4_raw(k)*100/subSum;
    CAP5(k) = CAP5_raw(k)*100/subSum;      CAP6(k) = CAP6_raw(k)*100/subSum;
end

tbl = table(sub', ses', subCat, CAP1_raw, CAP2_raw, CAP3_raw, CAP4_raw, CAP5_raw, CAP6_raw, CAP1, CAP2, CAP3, CAP4, CAP5, CAP6);
tbl.Properties.VariableNames = {'subject', 'session', 'group', 'raw CAP 1', 'raw CAP 2', 'raw CAP 3', 'raw CAP 4', 'raw CAP 5', 'raw CAP 6', 'CAP 1', 'CAP 2', 'CAP 3', 'CAP 4', 'CAP 5', 'CAP 6'};
writetable(tbl, fullfile(outPath, 'uCAPs_YODA.xlsx'),"FileType","spreadsheet")



