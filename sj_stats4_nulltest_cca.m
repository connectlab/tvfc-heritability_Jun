run sj_hmm_setting

if strcmp(answer{2}, '4')
    answer{3} = '16';
elseif strcmp(answer{2}, '6')
    answer{3} = '36';
end
%% Canonical Correlations
% Y = behavioral.
load(fullfile(dir_output, 'brain_behavior', ['cca_dynamic_brain_' answer{3} 'vars_pca_K' num2str(K) '_cogonly.mat']), 'cca_r', 'X', 'Y', 'cca_selected_mode')

iter = 10000;

for i = 1:iter
    X_perm = randperm(size(X, 1));
    [A,B,cca_r_perm(i, :),U,V,stats] = canoncorr(X(X_perm, :), Y); % X = behavioral (N x d1; N = subjects), Y = brain states (N x d2)
    cca_p_perm(i, :) = stats.pF;
    clear X_perm A B U V stats
end






