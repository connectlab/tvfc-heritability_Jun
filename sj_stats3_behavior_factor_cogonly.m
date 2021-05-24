%% Set up
clear; clc;

dir_scripts = pwd;
cd ..;
dir_home = pwd; cd (dir_scripts);
addpath(genpath(fullfile(dir_home, 'toolbox')));
dir_output = fullfile(dir_scripts, 'behavior'); mkdir(dir_output);

load twinid
subid_select = vertcat(pair_MZ(:), pair_DZ(:), pair_NT(:), pair_Unrelated(:));

%% Behavioral setting
load nptest_AgeAdj

num_idx = ismember(cell2mat(nptest_data(:, 1)), subid_select);
nptest_data = nptest_data(num_idx, :);

%%
bhv_domain = nptest_data(1, 5:end);
bhv_nptest = nptest_data(2, 5:end); bhv_nptest = replace(bhv_nptest, '_', '-');

% Classify the bhv data based on its domain
[domains domains_ind] = unique(bhv_domain, 'stable');
for i = 1:length(domains)
    eval(['bhv_data_' domains{i} '= bhv_data(:, ismember(bhv_domain, domains{i}));']);
end

domains = {'Cognition'};                    

% remove subjects with missing bhv data.
nanlist = sum(isnan(bhv_data_Cognition), 2);
for i = 1:length(nanlist)
    if nanlist(i)>0
        nanlists(i) = 1;
    else
        nanlists(i) = 0;
    end
end
nanlists = logical(nanlists);
bhv_data_Cognition = bhv_data_Cognition(~nanlists, :);

%% Data 
%centering each vector
for i = 1:size(bhv_data_Cognition, 2)
    bhv_data_centered(:, i) = (bhv_data_Cognition(:, i) - mean(bhv_data_Cognition(:, i)))/std(bhv_data_Cognition(:, i));
end

%% PCA of bhv

[PCA_COEFF, PCA_SCORE, PCA_LATENT, PCA_TSQUARED, PCA_EXPLAINED, PCA_MU] = pca(bhv_data_centered, 'Centered', false);
NumFactor = sum(PCA_LATENT>1);
VarExplained = sum(PCA_EXPLAINED(1:NumFactor));

for i = 1:length(bhv_nptest)
    pclabel{i} = ['PC' num2str(i)];
end

for i = 1:NumFactor
    PCA_influence{i, 1} = bhv_nptest{idx(1)}; 
    PCA_influence{i, 2} = PCA_COEFF(idx(1), i); 
    PCA_influence{i, 3} = bhv_nptest{idx(2)}; 
    PCA_influence{i, 4} = PCA_COEFF(idx(2), i); 
 
    clear M idx
end

save(fullfile(dir_output, 'bhv_pca_cogonly.mat'), 'PCA_*', 'NumFactor', 'VarExplained')

%% Factor analysis of bhv
[Factor_loading Factor_loadingvar Factor_T Factor_stats Factor_Score_regression] = factoran(bhv_data_centered, NumFactor, 'rotate', 'promax', 'scores', 'regression');
[Factor_loading Factor_loadingvar Factor_T Factor_stats Factor_Score_Bartlett] = factoran(bhv_data_centered, NumFactor, 'rotate', 'promax', 'scores', 'Bartlett');

for i = 1:NumFactor
    factorlabel{i} = ['Factor ' num2str(i)];
end

for i = 1:NumFactor
    Factor_influence{i, 1} = bhv_nptest{idx(1)};
    Factor_influence{i, 2} = Factor_loading(idx(1), i); 
    Factor_influence{i, 3} = bhv_nptest{idx(2)}; 
    Factor_influence{i, 4} = Factor_loading(idx(2), i); 
    
    clear M idx
end

save(fullfile(dir_output, 'bhv_factoran_cogonly.mat'), 'Factor_*')

%% Testing the score measures
Scores = horzcat(Factor_Score_Bartlett,Factor_Score_regression);
Scores_coeff = corrcoef(Scores);
for i = 1:size(Scores, 2)
    if i <= size(Scores, 2)/2
        score_label{i} = ['Factor' num2str(i) '-B'];
    else
        score_label{i} = ['Factor' num2str(i-size(Scores, 2)/2) '-R'];
    end
end


