run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

load twinid
twins = {'MZ','DZ','NT','Unrelated'};
subid_ind = 1:NSub;

%% One-way ANOVA of sibling status (TP & FO)
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'DMNCON', 'DMNDAN')

testing_metrics = {'DMNCON_norm', 'DMNDAN_norm'};

sibling_pair_DMNCON_norm = [];
sibling_pair_DMNDAN_norm = [];

% Redefine the tangent point: subtract surrogate mean from all subjects.
load(fullfile(dir_scripts, ['HMMrun_K' num2str(K) '_surrogate_mean.mat']))

DMNCON = DMNCON-DMNCON_surrogates_mean;
DMNDAN = DMNDAN-DMNDAN_surrogates_mean;


% The similarity (= Euclidean distance) between pairs.
sibling_pair_label = [];

for t = 1:length(twins)
    sibling_pair = eval(['pair_' twins{t} ';']);
    for s = 1:size(sibling_pair, 1)
        pair1 = subid_ind(ismember(subid, sibling_pair(s, 1)));
        pair2 = subid_ind(ismember(subid, sibling_pair(s, 2)));
        
        dmncon_pair1_nan = isnan(DMNCON(pair1, :)); dmncon_pair2_nan = isnan(DMNCON(pair2, :)); dmncon_nan = (dmncon_pair1_nan+dmncon_pair2_nan)>0;
        dmndan_pair1_nan = isnan(DMNDAN(pair1, :)); dmndan_pair2_nan = isnan(DMNDAN(pair2, :)); dmndan_nan = (dmndan_pair1_nan+dmndan_pair2_nan)>0;

        temp_DMNCON_norm(s, 1) = norm(DMNCON(pair1, ~dmncon_nan)-DMNCON(pair2, ~dmncon_nan));
        temp_DMNDAN_norm(s, 1) = norm(DMNDAN(pair1, ~dmndan_nan)-DMNDAN(pair2, ~dmndan_nan));
        
        clear pair1 pair2 *_nan
    end
    
    eval(['DMNCONnorm_' twins{t} '=temp_DMNCON_norm;']);
    eval(['DMNDANnorm_' twins{t} '=temp_DMNDAN_norm;']);
    
    sibling_pair_DMNCON_norm = vertcat(sibling_pair_DMNCON_norm, temp_DMNCON_norm);
    sibling_pair_DMNDAN_norm = vertcat(sibling_pair_DMNDAN_norm, temp_DMNDAN_norm);
      
    sibling_pair_label = vertcat(sibling_pair_label, repmat(twins(t), size(sibling_pair, 1), 1));
    clear sibling_pair temp_*
end

%% Run one-way ANOVA of sibling status
for m = 1:length(testing_metrics)
    [anova_p.(testing_metrics{m}),anova_tab.(testing_metrics{m}) ,anova_stats.(testing_metrics{m})] = anova1(eval(['sibling_pair_' testing_metrics{m}]), sibling_pair_label, 'off');
    title(['One-way ANOVA of sibling status: ' testing_metrics{m}])
    close(figure(1))
end


