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
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'Yeo_NewmanQ', 'DMNFPN')
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'FO_SubState', 'TP_Sub_ord_vector')

testing_metrics = {'Q','FC','FO','TP'};


% Redefine the tangent point: subtract surrogate mean from all subjects.
load(fullfile(dir_scripts, ['HMMrun_K' num2str(K) '_surrogate_mean.mat']))

FC = DMNFPN - FC_surrogates_mean; clear DMNFPN
Q = Yeo_NewmanQ-Q_surrogates_mean; clear Yeo_NewmanQ
FO = FO_SubState - FO_surrogates_mean; clear FO_SubState

ct = 0;
for i = 1:K
    for j = 1:K
        if i == j
        else ct = ct+1; 
            TP_surrogates_mean_vector(1, ct) = TP_surrogates_mean(i, j); 
        end
    end
end

for s = 1:size(TP_Sub_ord_vector, 1)
    TP(s, :) = TP_Sub_ord_vector(s, :) - TP_surrogates_mean_vector;
end
clear TP_Sub_ord_vector


% The similarity (= Euclidean distance) between pairs.
for m = 1:length(testing_metrics)
    sibling_label = [];
    state_label = [];
    
    eval(['sibling_pair_' testing_metrics{m} '= [];']);
    temp_metrics = eval(testing_metrics{m});
    
    for kk = 1:size(temp_metrics, 2)
        for g = 1:length(groups)
            sibling_pair = eval(['pair_' groups{g}]);
                      
            for s = 1:size(sibling_pair, 1)
                pair1 = find(ismember(subid, sibling_pair(s, 1))); 
                pair2 = find(ismember(subid, sibling_pair(s, 2))); 
                
                temp_norm_sibling(s, 1) = norm(temp_metrics(pair1, kk)- temp_metrics(pair2, kk));
                
                clear pair1 pair2
            end
            
            eval(['sibling_pair_' testing_metrics{m} ' = vertcat(sibling_pair_' testing_metrics{m} ', temp_norm_sibling);']);
            
            if kk == 1
                sibling_label = vertcat(sibling_label, repmat(groups(g), size(sibling_pair, 1), 1));
            end
            
            clear sibling_ppair s temp_norm_sibling
        end
        
        if strcmp(testing_metrics{m}, 'TP')
            state_label = vertcat(state_label, repmat({['TPcomp ' num2str(kk)]}, size(sibling_label, 1), 1));
        else
            state_label = vertcat(state_label, repmat({['K ' num2str(kk)]}, size(sibling_label, 1), 1));
        end
    end
    
    sibling_label = repmat(sibling_label, size(temp_metrics, 2), 1);
    
    eval(['state_label_' testing_metrics{m} '= state_label;']);
    eval(['sibling_label_' testing_metrics{m} '= sibling_label;']);
    
    clear sibling_label state_label temp_metrics
end

%% Two-way ANOVA analysis of sibling similarity x brain states.

for m = 1:length(testing_metrics)
    age_diff = vertcat(age_MZ(:, 1)-age_MZ(:, 2), age_DZ(:, 1)-age_DZ(:, 2), age_NT(:, 1)-age_NT(:, 2), age_Unrelated(:, 1)-age_Unrelated(:, 2));
    age_diff = repmat(age_diff, size(eval(testing_metrics{m}), 2), 1);
    data = eval(['sibling_pair_' testing_metrics{m}]);
    sibling_label = eval(['sibling_label_' testing_metrics{m}]);
    state_label = eval(['state_label_' testing_metrics{m}]);
    
    [anovan_p.(testing_metrics{m}),anovan_tab.(testing_metrics{m}),anovan_stats.(testing_metrics{m})] = anovan(data, {sibling_label, state_label, age_diff},'Continuous', 3, 'model','interaction','varnames',{'sibling status','State K', 'age_diff'});
     
    clear data sibling_label state_label age_diff 
end
