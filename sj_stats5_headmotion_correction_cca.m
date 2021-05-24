run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

if strcmp(answer{2}, '4')
    answer{3} = '16';
elseif strcmp(answer{2}, '6')
    answer{3} = '36';
end

%% Head motion and subjects 
load(fullfile(dir_home, 'scripts2_hmm_main', 'headmotion_FD.mat'), ['meanFD_subj_state_' num2str(K)])
load(fullfile(dir_scripts, 'behavior', 'bhv_data_cogonly.mat'), 'nanlists')

load twinid
subid_select = vertcat(pair_MZ(:), pair_DZ(:), pair_NT(:), pair_Unrelated(:));
subid_idx = ismember(subid, subid_select);
subid = subid(subid_idx);
subid = subid(~nanlists);

%%
FD = eval(['meanFD_subj_state_' num2str(K)]);
FD = nanmean(FD, 2);
FD = FD(subid_idx);
FD = FD(~nanlists);

load(fullfile(dir_output, 'brain_behavior', ['dynamic_brain_' answer{3} 'vars_pca_K' num2str(K) '.mat']), 'brain_data_heritable')
isnan_X = sum(isnan(brain_data_heritable), 2)>0;
FD = FD(~isnan_X);

%% Canonical Correlations
% Y = behavioral.
load(fullfile(dir_output, 'brain_behavior', ['cca_dynamic_brain_' answer{3} 'vars_pca_K' num2str(K) '_cogonly.mat']), 'cca_r', 'X', 'Y')

mdl_FD = fitlm(X, FD);                                                     % get residuals of X (brain) with head motion effect adjusted
R_FD = mdl_FD.Residuals(:, 4);                                             % get standardized FD model Residuals
R_FD = table2array(R_FD);

iter = 10000;

for i = 1:iter
    X_perm = randperm(size(X, 1));
    [A,B,cca_r_perm_FD(i, :),U,V,stats] = canoncorr(R_FD(X_perm, :), Y); % X = behavioral (N x d1; N = subjects), Y = brain states (N x d2)
    cca_p_perm_FD(i, :) = stats.pF;
    clear  A B U V stats

end



