%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Find best K ---
% mCAPs script based on the work of Farnaz Delavari (mCAP-main)
% testing range of K values to find best number of cluster
%
%
% Created by : Camille Serquet 
% Creation : 05.2025
% Last modification : 16.07.2025
% MATLAB version : R2022b
% SPM version : 12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*** START ***')

%% SETTING ENVIRONMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spmPath = fullfile('/srv/beegfs/scratch/users/s/serquet/spm12');            %[EDIT]: path to SPM toolbox
randomVecPath = '/srv/beegfs/scratch/users/s/serquet/randomVec.mat';        %[EDIT]: vector containing random order 
dataPath = fullfile('/srv/beegfs/scratch/users/s/serquet/dataV2.mat');      %[EDIT]: path to data.mat
outPath = fullfile('/srv/beegfs/scratch/users/s/serquet/outAllV2');         %[EDIT]: path to folder to store data

addpath(genpath(spmPath));                                                  %add path of SPM toolbox

if ~isfolder(outPath);  mkdir(outPath); end                                 %if output folder doesn't exist, create it

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kRange = (4:4);                                                             %[EDIT]: k range to be tested
maxIteration = 6;                                                           %[EDIT]: max number of iteration
convVal = 0.005;                                                            %[EDIT]: distance define as convergence across seeds
THR = 1;                                                                    %[EDIT]: threshold above which to select frames
Tmotion = 5;                                                                %[EDIT]: threshold for FD to scrub frames [t-1, t+1]
nRep = 100;                                                                 %[EDIT]: number of repetitions for kmeans
nIter = 300;                                                                %[EDIT]: number of iterations for kmeans
options = [];                                                               %[EDIT]: processing options for kmeans

SM = [1,0];                                                                 %sign matrix to retain frame - [1 0] activation - [0 1] deactivation
Pp = 100;                                                                   %percentage of positive-valued voxels to retain for clustering
Pn = 100;                                                                   %percentage of negative-valued voxels to retain for clustering
Pcc = 80;                                                                   %percentage of frames to use in each fold of consensus clustering

%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(dataPath);
disp(" >> Loading data : completed.")
disp(['  >> dimension of TC : ', num2str(size(TC))]);
disp(['  >> dimension of TC cell : ', num2str(size(TC{1,1}))]);
disp(['  >> dimension of FD : ', num2str(size(FD))]);
disp(['  >> dimension of seedmask : ', num2str(size(seedmask))]);
%% DEVIDE DATA INTO TRAIN AND TEST FRAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsub = size(TC,2);                                                          %get number of sessions 

%----------randomize order of sessions-------------------------------------
if exist(randomVecPath,'file')                                              %if a randomization already exists
    load(randomVecPath)                                                     %load the random vector
    if numel(randOrder) ~= nsub                                             %check if the size is the same as the current TC and if not:
        error("Saved random vector doesn't match size of current TC")       %display error and stop code
    end
else                                                                        %if a randomization doesn't already exist
    randOrder=randperm(nsub);                                               %create random vector of size nsub
    save(randomVecPath,"randOrder",'-mat')                                  %save the vector
end

%----------reorder TC and FD-----------------------------------------------
TCnew = TC(randOrder);
FDnew = FD(:, randOrder);
%-----------get half-set value---------------------------------------------
n = floor(nsub/2);
%-----------dividing TC into 2 subsets-------------------------------------
TCa = TCnew(1,1:n);
TCb = TCnew(1, n+1:nsub);
%-----------dividing FD into 2 subsets-------------------------------------
FDa = FDnew(:,1:n);
FDb = FDnew(:, n+1:nsub);
%-----------clear unnecessary variables------------------------------------
clear TC FD TCnew FDnew

%% INVESTIGATE K RANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kIt = 1:numel(kRange)
    %----------initialization of investigation-----------------------------
    ik = kRange(kIt);                                                       %k to be tested
    
    seedVol{1,1} = logical(seedmask);                                       %first seed volume = seedmask as single volume

    clear meanDis                                                           %clear variable to avoid errors
    
    isconverged = 0;                                                        %set initial convergence to 0
    it = 1;                                                                 %set iteration to 1

    [indAll, issAll, idx, Cm, SUMD, frameAssign, CAPnorm] = deal(cell(maxIteration, 1)); %initiate all variables
    [emptk, maxSUMD, meanDis] = deal(nan(maxIteration,1)); 
       

    disp(" >> Initialization : completed.")

    while isconverged==0 && it<maxIteration                                 %while iteration is not max and distance is not convergent
        disp(['>> Starting iteration n° ', num2str(it), '...']);

        signMatrix = repmat(SM, size(seedVol{it,1},2),1);                   %repete the sign matrix with number of seed volumes
        
        %----------finc activity in test and train frames------------------
        if it == 1                                                          %if we are doing the first iteration
            Tr = 0.5;                                                       %set the threshold to 0.5 
            [XonB,~,~,~]= CAPfindActivity(TCb,seedVol{it,1},Tr,FDb,Tmotion,signMatrix,seedmask, "CAP"); %find CAP activity
            [XonA,~,indAll{it,1},issAll{it,1}]= CAPfindActivity(TCa,seedVol{it,1},Tr,FDa,Tmotion,signMatrix,seedmask, "CAP");

        else                                                                %for all iterations after the first one
            [XonB,~,~,~]= CAPfindActivity(TCb,seedVol{it,1},THR,FDb,Tmotion,signMatrix,seedmask, "mCAP"); %find uCAP activity
            [XonA,~,indAll{it,1},issAll{it,1}]= CAPfindActivity(TCa,seedVol{it,1},THR,FDa,Tmotion,signMatrix,seedmask, "mCAP");
        end
        
        disp(['  >> dimension of XonB: ', num2str(size(XonB))]);
        disp(['  >> dimension of XonA: ', num2str(size(XonA))]);
        disp(['  >> dimension of indAll{it}: ', num2str(size(indAll{it}))]);
        disp(['  >> dimension of issAll{it}: ', num2str(size(issAll{it}))]);


        testFrames = cell2mat(XonB);                                        %test frames as matrix
        clear XonB;
        
        %----------k-means-------------------------------------------------
        [idx{it},Cm{it},SUMD{it},D{it}] = kmeans(cell2mat(XonA)',ik,'Options',options,'distance','correlation','replicates',nRep,'empty','drop','maxiter',nIter);
        trainFrames = cell2mat(XonA);                                       %train frames as matrix

        %----------test frames---------------------------------------------
        [emptk(it,1), maxSUMD(it,1), frameAssign{it},~] = tframe(testFrames,Cm{it});

        %----------normalize CAPs------------------------------------------
        meanC = zeros(size(trainFrames,1), ik);                             %initialize means matrix
        stdC = zeros(size(trainFrames,1), ik);                              %initialize std matrix
        for ii = 1:ik
            meanC(:,ii) = mean(trainFrames(:,idx{it}==ii),2,'omitnan');     %get mean
            stdC(:,ii) = std(trainFrames(:,idx{it}==ii),[],2,'omitnan');    %get std
        end


        CAPnorm{it} = meanC./stdC;                                          %normalize CAP value
        clear meanC stdC;
        
        %----------new seed division---------------------------------------
        newseed = CAPnorm{it}.*seedmask;                                    %divide seedmask into volumes
        it = it+1;                                                          %next iteration
        seedVol{it,1} = newseed;                                            %store seed volumes
        clear newseed;                                                      %clear variable
        
        %-----------assess convergence status------------------------------
        if it>2
            mnew = seedVol{it,1};                                           %current seed volumes
            mnew = mnew(seedmask==1,:);
            mold = seedVol{it-1,1};                                         %previous seed volumes
            mold = mold(seedmask==1,:);

            [finalDis] = compareSeeds(mnew, mold);                          %compare new with old seed
            meanDis(it,1) = mean(finalDis);                                 %calculate the mean distance btw seeds
                    
            if meanDis(it)<=convVal                                         %if the mean distance is lower than threshold
                isconverged = 1;                                            %mark it as convergent
            end
        end    
    end

%----------save variables to struct----------------------------------------
M.frameAssign = frameAssign;    M.CAPnorm = CAPnorm;   M.Cm = Cm;
M.D = D;                        M.idx = idx;           M.indAll = indAll;
M.issAll = issAll;              M.meanDis = meanDis;   M.seedVol = seedVol;
M.SUMD = SUMD;                  M.lastXON = XonA;
M.lastTestFrames = testFrames;  M.emptk = emptk;       M.maxSUMD = maxSUMD; 

meanMaxSumD = mean(maxSUMD, 'omitnan');

clear frameAssign CAPnorm Cm D idx indAll issAll seedVol SUMD 
clear XonA trainFrames testFrames maxSUMD

%----------saving output for i-th k----------------------------------------
savePath = fullfile(outPath,['TTmCAP_conv',char(num2str(isconverged)),'_K',char(string(ik))]);
mkdir(savePath);

saveName = fullfile(savePath,['mCAPTT_K',char(num2str(ik)),'.mat']);
save(saveName ,'emptk','meanMaxSumD','M','meanDis','it','-v7.3'); 

clear M it;

%----------write meanMaxSumD-----------------------------------------------
disp(['  >> meanMaxSumD : ', num2str(meanMaxSumD)])



end

disp('*** END ***')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xonp,p,ind,idxSeeds,XonPscrub,S] = CAPfindActivity(tc, seed,T,FDall,motTHR,SignMat, seedmask, CAPtype)
    %-----------initialization---------------------------------------------
    idxSeeds = nan(size(FDall,1),size(FDall,2),size(seed,2));               %initiate matrix storing active frames per seed
    flag = FDall>motTHR;                                                    %get frames exceeding motion threshold
    ind.scrubbed = logical(flag);                                           %save frames to exclude
    ind.kept.active = false(size(FDall,1), size(FDall,2));                  %initiate matrix of activ frames
    ind.scrubbedandactive = false(size(FDall,1), size(FDall,2));            %initiate matrix of scrubbed and activ frames
    flag = num2cell(flag, 1);                                               %transform to cell array
    pScrubbed = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag,'un',0));  %calculate percentage of scrubbed frames per subject
    
    %------------copmpute seed timecourse----------------------------------
    for i = 1:size(seed,2)
        %----------for CAP type (first it)---------------------------------
        if CAPtype == "CAP"
            seedA = logical(seed(:,i));                                     %get seed volume
            SignM = SignMat(i,:);                                           %get sign matrix for this seed
            if SignM(1)                                                     %if the sign is for activation - calculate seed timecourse
                S = cellfun(@(x) zscore(mean(x(:,seedA==1),2)), tc, 'UniformOutput', 0);        
            elseif SignM(2)                                                 %if the sign matrix is for deactivation - calculate seed timecourse
                S = cellfun(@(x) (-1)*zscore(mean(x(:,seedA==1),2)), tc, 'UniformOutput', 0);
            else                                                            %otherwise, error message
                errordlg('PROBLEM WITH SIGN MAT !!!');
            end

        %-----------for mCAP (after first it)------------------------------   
        elseif CAPtype == "mCAP"
            weightCAP = seed(:,i);                                          %get seed volume
  		    weightseed = weightCAP(seedmask==1,:);                          %extract seed voxels inside mask
   		    maxv = max(weightseed, [], 2);                                  %get max values in the masked seed

            maskforweight = zeros(size(weightseed));                        %initiate mask for weight seed
            for iv = 1:size(weightseed,1)                                   %for each voxel in weighted seed
                maskforweight(iv,:) = (weightseed(iv,:)==maxv(iv,1));       %keep only the column where weight is maximum
            end
  		    
  		    weightvalue = weightseed.*maskforweight;                        %keep only the max-weighted voxels assigne to the seed
           
            powerofsum1 = 1;                                                %power controlling sharpness of voxel weights on signal
            powerofsum2 = 1;                                                %power controlling the normalization factor computation
            weightpower = weightvalue.^powerofsum1;                         %adjust infulence of weight
            denim = sum((abs(weightvalue)).^powerofsum2);                   %normalize final value
            
            SignM = SignMat(i,:);                                           %get sign matrix for this seed
            if SignM(1)                                                     %if the sign is for activation - calculate seed timecourse
                S = cellfun(@(x) zscore((sum((x(:,seedmask).*weightpower'),2))./denim), tc, 'UniformOutput', 0);
            elseif SignM(2)                                                 %if the sign is for deactivation - calculate seed timecourse
                S = cellfun(@(x) (-1)*zscore((sum((x(:,seedmask).*weightpower'),2))./denim), tc, 'UniformOutput', 0);
            else                                                            %otherwise, error
                errordlg('PROBLEM WITH SIGN MAT !!!');
            end
        else                                                                %if CAP type is wrong, error
            error("NOPE")
        end
        
        %-----------activation threshold and scrubbing---------------------
        xindp = cellfun(@(x) x>T, S, 'un', 0);                                              %find timepoints where signal exceeds threshold
        flagActive = cellfun(@(x,y) x & y,xindp, flag,'un',0);                              %find timepoints that are active ans exceed threshold
        pScrubActive = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flagActive,'un',0));     %get percentage of active and scrubbed timepoints
        xindp = cellfun(@(x,y) x & ~y, xindp,flagActive,'un',0);                            %get timepoints active and motion-free
        ind.scrubbedandactive = (ind.scrubbedandactive) | cell2mat(flagActive);             %store scrubbed and active frames
        ind.kept.active = (ind.kept.active) | cell2mat(xindp);                              %store active frames
        idxSeeds(:,:,i) = cell2mat(xindp);                                                  %update per-seed tracking of active frames  
    end
    
    %----------extract data from active frames-----------------------------
    Xonp = cellfun(@(x,y) x(y,:)', tc, num2cell(ind.kept.active,1), 'UniformOutput', 0);                        %keep active frames
    XonPscrub = cellfun(@(x,y) x(y,:)', tc, num2cell(logical(ind.scrubbedandactive), 1), 'UniformOutput', 0);   %keep active and scrubbed frames
    pActive = cell2mat(cellfun(@(x) size(x,2)/size(FDall,1)*100, Xonp, 'UniformOutput', 0));                    %compute percentage of active frames
    p = [pScrubbed; pScrubActive; pActive];                                                                     %combine motion scrubbed, scrubbed and active, active and kept frames
end


function [finalDis] = compareSeeds(seed1,seed2)
    seed1 = double(seed1);                                                  %seed 1
    seed2 = double(seed2);                                                  %seed 2

    DIS = pdist2(seed1',seed2','cosine');                                   %compute pairwise cosine distance btw seed 1 and seed 2

    IDX = munkres(DIS);                                                     %find optimal matching
    
    finalDis = zeros(size(DIS,1), 1);                                       %initiate final distance vector
    idxk = 1;                                                               %initiate iteration 
    for k = 1:size(DIS,1)                                                   %for all distances
        if sum(IDX(k,:)) == 1                                               %if an optimal match is found
            finalDis(idxk) = DIS(k,IDX(k,:));                               %the distance is this optimisation one
            idxk = idxk + 1;                                                %go to next iteration
        end
    end
end

    
function [assignment, cost] = munkres(costMat)

% costMat :     matrix where the cost of assigning job j to worker i
% assignment :  binary matrix
% cost :        total cost of the optimal assignment

    %----------- initialization--------------------------------------------
    assignment = false(size(costMat));                                      %initiate assign matrix with false values
    cost = 0;                                                               %initialise total assignment cost
    costMat(costMat~=costMat)=Inf;                                          %replace NaN with Inf (insure proper computation)
    validMat = costMat<Inf;                                                 %determine which elements are valid (not Inf)
    validCol = any(validMat);                                               %identity valid column
    validRow = any(validMat,2);                                             %identify valid rows
    nRows = sum(validRow);                                                  %count nb of valid rows
    nCols = sum(validCol);                                                  %count nb of valid column
    n = max(nRows,nCols);                                                   %determine larger dimension
    if ~n                                                                   %if there are no valid assignement (all values are Inf)
        return                                                              %return early
    end
        
    dMat = zeros(n);                                                        %initiate a square zero matrix
    dMat(1:nRows,1:nCols) = costMat(validRow,validCol);                     %fill the first row/col with valid cost matrix values
    
    %----------substract row minimum---------------------------------------
    
    dMat = bsxfun(@minus, dMat, min(dMat,[],2));                            %substract row minimum from each row
    
    %----------star (*) zeros----------------------------------------------
    zP = ~dMat;                                                             %get zeros in the matrix
    starZ = false(n);                                                       %create false matrix to store *0
    
    while any(zP(:))                                                        %while there are still unprocessed 0
        [r,c]=find(zP,1);                                                   %find the first available 0
        starZ(r,c)=true;                                                    %star the 0 ___ 0 -> *0
        zP(r,:)=false;                                                      %remove 0s in the same row
        zP(:,c)=false;                                                      %remove 0s in the same column
    end
    
    %----------cover *0----------------------------------------------------
    
    while 1                                                                 %while it's the case - until we stop the loop
        primeZ = false(n);                                                  %initiate '0
        cColumn = any(starZ);                                               %get column containing a *0
        if ~any(~cColumn)                                                   %if all columns are covered
            break                                                           %break while loop
        end

        cRow = false(n,1);                                                  %initiate covered row array
        while 1                                                             %while it's the case - until we stop the loop
           zP(:) = false;                                                   %put all values to false
            zP(~cRow,~cColumn) = ~dMat(~cRow,~cColumn);                     %find uncovered zeros
            step = 6;                                                       %set up step to 6
            while any(any(zP(~cRow,~cColumn)))                              %while there are uncovered 0
                [uZr,uZc] = find(zP,1);                                     %find the first uncovered 0
                primeZ(uZr,uZc) = true;                                     %prime it ___ 0 -> '0
                stz = starZ(uZr,:);                                         %find *0 in the row
                if ~any(stz)                                                %if there's no *0 in this row
                    step = 5;                                               %set up step to 5
                    break;                                                  %break the loop
                end
                cRow(uZr) = true;                                           %cover the row of the first uncovered 0
                cColumn(stz) = false;                                       %uncover the column containing the *0
                zP(uZr,:) = false;                                          %remove processed 0
                zP(~cRow,stz) = ~dMat(~cRow,stz);                           %rewrite matrix with processed values
            end
            
            if step == 6                                                    %if the step is set to 6
                M = dMat(~cRow,~cColumn);                                   %get uncovered matrix
                minval = min(min(M));                                       %find minimum value
                if minval == inf                                            %if the min value is infiny
                    return                                                  %stop loop
                end

                dMat(cRow,cColumn) = dMat(cRow,cColumn) + minval;           %add min val to covered row
                dMat(~cRow,~cColumn) = M - minval;                          %substract min val to uncovered row
            
            else                                                            %if we don't need step 6
                break                                                       %break loop
            end
        end
        
        rowZ1 = starZ(:, uZc);                                              %find *0 in the column of uncovered '0
        starZ(uZr, uZc) = true;                                             %star the uncovered '0

        while any(rowZ1)                                                    %while we still have *0
            starZ(rowZ1, uZc) = false;                                      %unstar it
            uZc = primeZ(rowZ1, :);                                         %find the '0 in its row
            uZr = rowZ1;                                                    %new column containing *0
            rowZ1 = starZ(:, uZc);                                          %find the *0 in the new column
            starZ(uZr, uZc) = true;                                         %star the '0
        end
    end

    %----------cost of assignment------------------------------------------
    assignment(validRow,validCol) = starZ(1:nRows, 1:nCols);                %extract optimal assignment from the *0s
    cost = sum(costMat(assignment));                                        %compute total assignment cost
end


function [emptK,maxSumD,fAssign,sumD] = tframe(frame,Cma)
%assign a set of frames to the closest centroide/clusteer
    
    %----------initialization of variables---------------------------------
    k = size(Cma,1);                                                        %initiate number of centroids/clusters
    f = size(frame,2);                                                      %initiate number of frames
    fAssign = zeros(f,1);                                                   %initiate cluster assignment for each frame
    disFassign = zeros(f,1);                                                %initiate distance to the assigned cluster
    
    %----------assign frames to closest centroid---------------------------
    for im = 1:f                                                            %for each frames
        disMat = zeros(k,1);                                                %initiate a new distance matrix
       for ic = 1:k                                                         %for all centroids
            disMat(ic,1) = pdist2(frame(:,im)',Cma(ic,:),'correlation');    %assign each frame to the closest centroid       
       end
       
       targetcap = find(disMat==min(disMat));                               %find closest centroid
       fAssign(im,1) = targetcap(1);                                        %assign to closest
       disFassign(im,1) = min(disMat);                                      %store the distance   
    end

    %----------get sum of distances per clusters---------------------------
    sumD = zeros(k,1);                                                      %initiate vector to store sum of distance    
    for ic = 1:k                                                            %for each centroids
        sumD(ic,1)=sum(disFassign(fAssign==ic));                            %sum distances per clusters
    end

    %----------get empty cluster and max distance--------------------------
    emptK = k-size(unique(fAssign),1);                                      %get nb of empty cluster
    maxSumD = max(sumD);                                                    %cluster with max of total distance
end
















