run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

load twinid
twins = {'MZ','DZ','NT','Unrelated'};
subid_ind = 1:NSub;

%% One-way ANOVA of sibling status (TP & FO)
load(fullfile(dir_scripts, 'headmotion_FD.mat'), 'meanFD_subj_state*')
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'TP_Sub_ord', 'FO_SubState')
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'Yeo_NewmanQ', 'DMNFPN')
FC = DMNFPN; clear DMNFPN
FD = eval(['meanFD_subj_state_' num2str(K)]);

testing_metrics = {'TP_norm', 'FO_norm', 'Q_norm', 'FC_norm', 'FD_norm'};

% TP & FO corr
twin_pair_TP_norm = [];
twin_pair_FO_norm = [];
twin_pair_Q_norm = [];
twin_pair_FC_norm = [];
twin_pair_FD_norm = [];

twin_pair_label = [];

for t = 1:length(twins)
    twin_pair = eval(['pair_' twins{t} ';']);
    for s = 1:size(twin_pair, 1)
        pair1 = subid_ind(ismember(subid, twin_pair(s, 1)));
        pair2 = subid_ind(ismember(subid, twin_pair(s, 2)));
        
        fo_pair1_nan = isnan(FO_SubState(pair1, :)); fo_pair2_nan = isnan(FO_SubState(pair2, :)); fo_nan = (fo_pair1_nan+fo_pair2_nan)>0;
        q_pair1_nan = isnan(Yeo_NewmanQ(pair1, :)); q_pair2_nan = isnan(Yeo_NewmanQ(pair2, :)); q_nan = (q_pair1_nan+q_pair2_nan)>0;
        fc_pair1_nan = isnan(FC(pair1, :)); fc_pair2_nan = isnan(FC(pair2, :)); fc_nan = (fc_pair1_nan+fc_pair2_nan)>0;
        fd_pair1_nan = isnan(FD(pair1, :)); fd_pair2_nan = isnan(FD(pair2, :)); fd_nan = (fd_pair1_nan+fd_pair2_nan)>0;

        temp_TP_norm(s, 1) = norm((TP_Sub_ord{pair1}- TP_Sub_ord{pair2}), 2);
        temp_FO_norm(s, 1) = norm(FO_SubState(pair1, ~fo_nan)-FO_SubState(pair2, ~fo_nan));
        temp_Q_norm(s, 1) = norm(Yeo_NewmanQ(pair1, ~q_nan)-Yeo_NewmanQ(pair2, ~q_nan));
        temp_FC_norm(s, 1) = norm(FC(pair1, ~fc_nan)-FC(pair2, ~fc_nan));
        temp_FD_norm(s, 1) = norm(FD(pair1, ~fd_nan)-FC(pair2, ~fd_nan));
        
        clear pair1 pair2 *_nan
    end
    
    eval(['TPnorm_' twins{t} '=temp_TP_norm;']);
    eval(['FOnorm_' twins{t} '=temp_FO_norm;']);
    eval(['Qnorm_' twins{t} '=temp_Q_norm;']);
    eval(['FCnorm_' twins{t} '=temp_FC_norm;']);
    eval(['FDnorm_' twins{t} '=temp_FD_norm;']);
    
    
    twin_pair_TP_norm = vertcat(twin_pair_TP_norm, temp_TP_norm);
    twin_pair_FO_norm = vertcat(twin_pair_FO_norm, temp_FO_norm);
    twin_pair_Q_norm = vertcat(twin_pair_Q_norm, temp_Q_norm);
    twin_pair_FC_norm = vertcat(twin_pair_FC_norm, temp_FC_norm);
    twin_pair_FD_norm = vertcat(twin_pair_FD_norm, temp_FD_norm);
    
    twin_pair_label = vertcat(twin_pair_label, repmat(twins(t), size(twin_pair, 1), 1));
    clear twin_pair temp_*
end


