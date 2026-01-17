clear; clc; close all;

%% Mask (Amygdala indices)
V = spm_vol('Mask file');
mask = spm_read_vols(V);
Amyg = find(mask == 1);

%% Load data
load('fMRI data sort by valence'); % group, SortData, ValenceOrder
load('Sorted valence') 

valence = ValenceOrder;

%% Build data matrix: [nAmygVox x (60*nSub)]
nSub = 20;
nPer = 60;

tmp = SortData(Amyg,:,1:nSub);                 % [nAmygVox x nPer x nSub]
DataA_sort_C = reshape(cat(2, tmp), numel(Amyg), nPer*nSub);

% Expand labels

valence_M = repmat(valence(:), nSub, 1);

%% Remove voxels (rows) that contain any NaN across samples
validRow = all(~isnan(DataA_sort_C), 2);
Amyg = Amyg(validRow);
X = zscore(DataA_sort_C(validRow,:), 0, 2);   % zscore within voxel across samples
Y = valence_M;

%% Univariate 
% X = mean(X, 1);         % Uncomment this line when perform univariate regression (mean across voxels)                       % 1 x (60*nSub)

%% Test regression performance with cross-validation (leave-1-subject-out)
kfold  = nSub;
lambda = 1e-3;

R_train = zeros(kfold,1);
R_test  = zeros(kfold,1);



for k = 1:kfold
    idxTest = false(1, nPer*nSub);
    idxTest((k-1)*nPer + (1:nPer)) = true;
    idxTrain = ~idxTest;

    Mdl = fitrlinear(X(:,idxTrain), Y(idxTrain), ...
        'Solver','sparsa','Regularization','lasso','Lambda',lambda, ...
        'ObservationsIn','columns','Learner','leastsquares');

    yTrainHat = predict(Mdl, X(:,idxTrain), 'ObservationsIn','columns');
    R_train(k) = corr(yTrainHat(:), Y(idxTrain));

    yHat = predict(Mdl, X(:,idxTest), 'ObservationsIn','columns');
    [R_test(k), P_test(k)] = corr(zscore(yHat(:)), zscore(Y(idxTest)));
    MSE(k) = immse(zscore(yHat(:)), zscore(Y(idxTest)));

    BetaMDL2(:,k)=Mdl.Beta'*cov(X(:,idxTrain)'); % Weight map
end

R_train_mean = mean(R_train);
R_test_mean  = mean(R_test);



