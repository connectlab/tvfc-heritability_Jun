run sj_hmm_setting
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

load twinid
twins = {'MZ','DZ','NT','Unrelated'};
subid_ind = 1:NSub;

subid_MZ = pair_MZ(:);
subid_DZ = pair_DZ(:);
subid_NT = pair_NT(:);
subid_Unrelated = pair_Unrelated(:);

%% Two-way ANOVA analysis of sibling status x brain states.
load(fullfile(dir_home, 'scripts2_hmm_main', 'headmotion_FD.mat'), 'meanFD_subj_state*')
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'Yeo_NewmanQ', 'DMNFPN')
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'FO_SubState', 'TP_Sub_ord_vector')

Q = Yeo_NewmanQ; clear Yeo_NewmanQ
FC = DMNFPN; clear DMNFPN
FO = FO_SubState; clear FO_SubState
TP = TP_Sub_ord_vector;
testing_metrics = {'Q','FC','FO','TP'};

for m = 1:length(testing_metrics)
    twin_pair_label_temp = [];
    state_label_indiv = [];
    
    eval(['twin_pair_' testing_metrics{m} '_indiv= [];']);
    temp_metrics = eval(testing_metrics{m});
    
    for i = 1:size(temp_metrics, 2)
        for t = 1:length(twins)
            twin_pair = eval(['pair_' twins{t}]);
            pair1 = subid_ind(ismember(subid, twin_pair(:, 1)));
            pair2 = subid_ind(ismember(subid, twin_pair(:, 2)));
            
            for s = 1:size(twin_pair, 1)
                temp_norm_twin(s, 1) = norm(temp_metrics(pair1(s), i)- temp_metrics(pair2(s), i));
            end
            
            eval(['twin_pair_' testing_metrics{m} '_indiv = vertcat(twin_pair_' testing_metrics{m} '_indiv, temp_norm_twin);']);
            
            if i == 1
                twin_pair_label_temp = vertcat(twin_pair_label_temp, repmat(twins(t), size(twin_pair, 1), 1));
            end
            clear temp_norm_twin twin_pair pair1 pair2
        end
        
        if strcmp(testing_metrics{m}, 'TP')
            state_label_indiv = vertcat(state_label_indiv, repmat({['TPcomp ' num2str(i)]}, size(twin_pair_label_temp, 1), 1));
        else
            state_label_indiv = vertcat(state_label_indiv, repmat({['K ' num2str(i)]}, size(twin_pair_label_temp, 1), 1));
        end
    end
    
    twin_pair_label_temp = repmat(twin_pair_label_temp, size(temp_metrics, 2), 1);
    
    eval(['state_' testing_metrics{m} '_indiv_label= state_label_indiv;']);
    eval(['twin_pair_' testing_metrics{m} '_indiv_label= twin_pair_label_temp;']);
    
    clear temp_metrics twin_pair_label_temp state_label_indiv
end

%%
for m = 1:length(testing_metrics)
    clear data sibling_label state_label
    data = eval(['twin_pair_' testing_metrics{m} '_indiv']);
    sibling_label = eval(['twin_pair_' testing_metrics{m} '_indiv_label']);
    state_label = eval(['state_' testing_metrics{m} '_indiv_label']);
    
    [anova2_p.(testing_metrics{m}),anova2_tab.(testing_metrics{m}),anova2_stats.(testing_metrics{m})] = anovan(data, {sibling_label, state_label},'model','interaction','varnames',{'sibling status','State K'}, 'display', 'off');
    multcompare(anova2_stats.(testing_metrics{m}));
        
    clear temp_* data sibling_label state_label
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova2_siblingsimilairity_x_state.mat']), 'anova2_*', 'twin_*', '*_label')


