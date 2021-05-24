run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

dir_bhv = fullfile(dir_scripts, 'behavior');
load(fullfile(dir_bhv, 'bhv_data_cogonly.mat'), 'nanlists')

load twinid
subid_select = vertcat(pair_MZ(:), pair_DZ(:), pair_NT(:), pair_Unrelated(:));
subid_idx = ismember(subid, subid_select);
subid = subid(subid_idx);
subid = subid(~nanlists);

mkdir(dir_output, 'brain_behavior')

%% brain Setting (dynamic brain measuers)
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'FO_SubState', 'TP_Sub_ord_vector', 'TP_coord')

for j = 1:K
    label_fo{j} = ['FO' num2str(j)];
end
for i = 1:size(TP_coord)
    label_tp{i} = ['TP' num2str(TP_coord(i,1)) num2str(TP_coord(i,2))];
end

label_heritable = horzcat(label_fo, label_tp);
brain_data_heritable = horzcat(FO_SubState, TP_Sub_ord_vector);

brain_data_heritable = brain_data_heritable(subid_idx, :);
brain_data_heritable = brain_data_heritable(~nanlists, :);

NSub = length(subid);
brain_data_heritable = brain_data_heritable-mean(brain_data_heritable, 'omitnan'); %centering
for i = 1:size(brain_data_heritable, 2)
    brain_data_heritable(:, i) = brain_data_heritable(:, i)/std(brain_data_heritable(:, i), 'omitnan'); %scaling
end

%% PCA of brain

[PCA_COEFF, PCA_SCORE, PCA_LATENT, PCA_TSQUARED, PCA_EXPLAINED, PCA_MU] = pca(brain_data_heritable, 'Centered', false);
NumFactor = sum(PCA_LATENT>1);
VarExplained = sum(PCA_EXPLAINED(1:NumFactor));

for i = 1:length(label_heritable)
    pclabel{i} = ['PC' num2str(i)];
end
    
for i = 1:NumFactor
    PCA_influence{i, 1} = label_heritable{idx(1)};
    PCA_influence{i, 2} = PCA_COEFF(idx(1), i);
    PCA_influence{i, 3} = label_heritable{idx(2)};
    PCA_influence{i, 4} = PCA_COEFF(idx(2), i);
    
    clear M idx
end

save(fullfile(dir_output, 'brain_behavior', ['dynamic_brain_' num2str(length(label_heritable)) 'vars_pca_K' num2str(K) '.mat']), 'PCA_*', 'NumFactor', 'pclabel', 'VarExplained', 'brain_data_heritable', 'label_heritable')




