run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

load twinid
twins = {'MZ','DZ','NT','Unrelated'};

if strcmp(answer{2}, '4')
    answer{3} = '16';
elseif strcmp(answer{2}, '6')
    answer{3} = '36';
end
%% Behavioral data Setting
dir_bhv = fullfile(dir_scripts, 'behavior');
load(fullfile(dir_bhv, 'bhv_data_cogonly.mat'), 'nanlists')
load(fullfile(dir_bhv, 'bhv_factoran_cogonly.mat'), 'Factor_Score_regression*') %Factor_Score_Bartlett Factor_Score_regression
load(fullfile(dir_bhv, 'bhv_pca_cogonly.mat'), 'NumFactor')

for i = 1:NumFactor
    label_factor{i} = ['F' num2str(i)];
end

subid = subid(~nanlists);
NSub = length(subid);

Factor_Score_regression = Factor_Score_regression(:, 1:NumFactor);

x_table = table(Factor_Score_regression, 'RowNames',cellstr(num2str(subid)));
x_table = splitvars(x_table, 'Factor_Score_regression', 'NewVariableNames', label_factor);

%% brain Setting (dynamic brain measuers)

load(fullfile(dir_output, 'brain_behavior', ['dynamic_brain_' answer{3} 'vars_pca_K' num2str(K) '.mat']))
brain_data_heritable = PCA_SCORE(:, 1:NumFactor);
label_heritable = pclabel;

%% Canonical Correlations
% Y = behavioral.
isnan_X = sum(isnan(brain_data_heritable), 2)>0;
X = brain_data_heritable(~isnan_X, :);
Y = Factor_Score_regression(~isnan_X, :);
[cca_A,cca_B,cca_r,cca_U,cca_V,cca_stats] = canoncorr(X, Y); % X = behavioral (N x d1; N = subjects), Y = brain states (N x d2)

cca_selected_mode = cca_stats.pF <.05;
cca_U_selected = cca_U(:, cca_selected_mode);
[cca_loading_bhv cca_loading_bhv_p]= corr(cca_U_selected, Factor_Score_regression(~isnan_X, :));
cca_V_selected = cca_V(:, cca_selected_mode);
[cca_loading_brain cca_loading_brain_p]= corr(cca_V_selected, brain_data_heritable(~isnan_X, :));

nmode = sum(cca_selected_mode);

save(fullfile(dir_output,'brain_behavior', ['cca_dynamic_brain_' answer{3} 'vars_pca_K' num2str(K) '_cogonly.mat']), 'cca*', 'X', 'Y')








