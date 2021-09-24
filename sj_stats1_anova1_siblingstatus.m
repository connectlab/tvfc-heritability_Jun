% Assuming that the tangent point set to compute the Euclidean distance
% bewteen the subject pairs would be n dimensional zeroes (e.g., 0 0 0 0),
% we extracted the surrogate mean obtained from the 50 simulated datasets
% from every individual, in order to make the statistics and tangent point
% comparable to the other version of statistical test.

run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

load twinid

%% One-way ANOVA of sibling status

load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'TP_Sub_ord', 'FO_SubState')
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'Yeo_NewmanQ', 'DMNFPN', 'DMNCON', 'DMNDAN')
FC = DMNFPN; clear DMNFPN
Q = Yeo_NewmanQ;

load(fullfile(dir_home, 'scripts2_hmm_main', 'headmotion_FD.mat'), 'FD_sub_K*') %not resampled
FD = eval(['FD_sub_K' num2str(K)]);

testing_metrics = {'TP_norm', 'FO_norm', 'Q_norm', 'FC_norm', 'FD_norm', 'DMNCON_norm', 'DMNDAN_norm'};


% Redefine the tangent point: subtract surrogate mean from all subjects.
dir_simhmmmar_stats = replace(dir_scripts, 'scripts2_hmm_main_stats', 'scripts5_hmm_simhmmmar_stats');
load(fullfile(dir_simhmmmar_stats, ['HMMrun_K' num2str(K) '_surrogate_mean.mat']))

for s = 1:length(TP_Sub_ord)
    TP_Sub_ord{s} = TP_Sub_ord{s} - TP_surrogates_mean;
end

FO_SubState = FO_SubState - FO_surrogates_mean;
FC = FC - FC_surrogates_mean;
Q = Q-Q_surrogates_mean;
DMNCON = DMNCON-DMNCON_surrogates_mean;
DMNDAN = DMNDAN-DMNDAN_surrogates_mean;


% The similarity (= Euclidean distance) between pairs.
sibling_pair_TP_norm = [];
sibling_pair_FO_norm = [];
sibling_pair_Q_norm = [];
sibling_pair_FC_norm = [];
sibling_pair_FD_norm = [];
sibling_pair_DMNCON_norm = [];
sibling_pair_DMNDAN_norm = [];

sibling_pair_label = [];

for t = 1:length(groups)
    sibling_pair = eval(['pair_' groups{t} ';']);
    
    for s = 1:size(sibling_pair, 1)
        pair1 = find(ismember(subid, sibling_pair(s, 1)));
        pair2 = find(ismember(subid, sibling_pair(s, 2)));
        
        fo_pair1_nan = isnan(FO_SubState(pair1, :)); fo_pair2_nan = isnan(FO_SubState(pair2, :)); fo_nan = (fo_pair1_nan+fo_pair2_nan)>0;
        q_pair1_nan = isnan(Q(pair1, :)); q_pair2_nan = isnan(Q(pair2, :)); q_nan = (q_pair1_nan+q_pair2_nan)>0;
        fc_pair1_nan = isnan(FC(pair1, :)); fc_pair2_nan = isnan(FC(pair2, :)); fc_nan = (fc_pair1_nan+fc_pair2_nan)>0;
        fd_pair1_nan = isnan(FD(pair1, :)); fd_pair2_nan = isnan(FD(pair2, :)); fd_nan = (fd_pair1_nan+fd_pair2_nan)>0;
        dmncon_pair1_nan = isnan(DMNCON(pair1, :)); dmncon_pair2_nan = isnan(DMNCON(pair2, :)); dmncon_nan = (dmncon_pair1_nan+dmncon_pair2_nan)>0;
        dmndan_pair1_nan = isnan(DMNDAN(pair1, :)); dmndan_pair2_nan = isnan(DMNDAN(pair2, :)); dmndan_nan = (dmndan_pair1_nan+dmndan_pair2_nan)>0;
       
        temp_TP_norm(s, 1) = norm((TP_Sub_ord{pair1}- TP_Sub_ord{pair2}), 2);
        temp_FO_norm(s, 1) = norm(FO_SubState(pair1, ~fo_nan)-FO_SubState(pair2, ~fo_nan));
        temp_Q_norm(s, 1) = norm(Q(pair1, ~q_nan)-Q(pair2, ~q_nan));
        temp_FC_norm(s, 1) = norm(FC(pair1, ~fc_nan)-FC(pair2, ~fc_nan));
        temp_FD_norm(s, 1) = norm(FD(pair1, ~fd_nan)-FD(pair2, ~fd_nan));
        temp_DMNCON_norm(s, 1) = norm(DMNCON(pair1, ~dmncon_nan)-DMNCON(pair2, ~dmncon_nan));
        temp_DMNDAN_norm(s, 1) = norm(DMNDAN(pair1, ~dmndan_nan)-DMNDAN(pair2, ~dmndan_nan));
        
        clear pair1 pair2 *_nan
    end
    
    eval(['TPnorm_' groups{t} '=temp_TP_norm;']);
    eval(['FOnorm_' groups{t} '=temp_FO_norm;']);
    eval(['Qnorm_' groups{t} '=temp_Q_norm;']);
    eval(['FCnorm_' groups{t} '=temp_FC_norm;']);
    eval(['FDnorm_' groups{t} '=temp_FD_norm;']);
    eval(['DMNCONnorm_' groups{t} '=temp_DMNCON_norm;']);
    eval(['DMNDANnorm_' groups{t} '=temp_DMNDAN_norm;']);
    
    sibling_pair_TP_norm = vertcat(sibling_pair_TP_norm, temp_TP_norm);
    sibling_pair_FO_norm = vertcat(sibling_pair_FO_norm, temp_FO_norm);
    sibling_pair_Q_norm = vertcat(sibling_pair_Q_norm, temp_Q_norm);
    sibling_pair_FC_norm = vertcat(sibling_pair_FC_norm, temp_FC_norm);
    sibling_pair_FD_norm = vertcat(sibling_pair_FD_norm, temp_FD_norm);
    sibling_pair_DMNCON_norm = vertcat(sibling_pair_DMNCON_norm, temp_DMNCON_norm);
    sibling_pair_DMNDAN_norm = vertcat(sibling_pair_DMNDAN_norm, temp_DMNDAN_norm);
    
    sibling_pair_label = vertcat(sibling_pair_label, repmat(groups(t), size(sibling_pair, 1), 1));
    clear sibling_pair temp_*
end


%% Visual Inspection of the main variables (TP, FO, FC, Q)

figure(1)
clf
set(figure(1), 'color', 'white', 'defaultAxesFontSize', 15)
tiledlayout(2, 2)

nexttile
Y = [mean(FOnorm_MZ, 'omitnan'), mean(FOnorm_DZ, 'omitnan'), mean(FOnorm_NT, 'omitnan'), mean(FOnorm_Unrelated, 'omitnan')];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [std(FOnorm_MZ, 'omitnan')/sqrt(length(FOnorm_MZ)), std(FOnorm_DZ, 'omitnan')/sqrt(length(FOnorm_DZ)), std(FOnorm_NT, 'omitnan')/sqrt(length(FOnorm_NT)), std(FOnorm_Unrelated, 'omitnan')/sqrt(length(FOnorm_Unrelated))];
er = errorbar(X,Y,Y_se/2,Y_se/2);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Fractional Occupancy')
ylabel('Euclidean distance')
axis square
hold off

nexttile
Y = [mean(TPnorm_MZ, 'omitnan'), mean(TPnorm_DZ, 'omitnan'), mean(TPnorm_NT, 'omitnan'), mean(TPnorm_Unrelated, 'omitnan')];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [std(TPnorm_MZ, 'omitnan')/sqrt(length(TPnorm_MZ)), std(TPnorm_DZ, 'omitnan')/sqrt(length(TPnorm_DZ)), std(TPnorm_NT, 'omitnan')/sqrt(length(TPnorm_NT)), std(TPnorm_Unrelated, 'omitnan')/sqrt(length(TPnorm_Unrelated))];
er = errorbar(X,Y,Y_se/2,Y_se/2);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Transition Probability')
ylabel('Euclidean distance')
axis square
hold off

nexttile
Y = [mean(FCnorm_MZ, 'omitnan'), mean(FCnorm_DZ, 'omitnan'), mean(FCnorm_NT, 'omitnan'), mean(FCnorm_Unrelated, 'omitnan')];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [std(FCnorm_MZ, 'omitnan')/sqrt(length(FCnorm_MZ)), std(FCnorm_DZ, 'omitnan')/sqrt(length(FCnorm_DZ)), std(FCnorm_NT, 'omitnan')/sqrt(length(FCnorm_NT)), std(FCnorm_Unrelated, 'omitnan')/sqrt(length(FCnorm_Unrelated))];
er = errorbar(X,Y,Y_se/2,Y_se/2);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('TV-FC DMN-FPN')
ylabel('Euclidean distance')
axis square
hold off

nexttile
Y = [mean(Qnorm_MZ, 'omitnan'), mean(Qnorm_DZ, 'omitnan'), mean(Qnorm_NT, 'omitnan'), mean(Qnorm_Unrelated, 'omitnan')];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [std(Qnorm_MZ, 'omitnan')/sqrt(length(Qnorm_MZ)), std(Qnorm_DZ, 'omitnan')/sqrt(length(Qnorm_DZ)), std(Qnorm_NT, 'omitnan')/sqrt(length(Qnorm_NT)), std(Qnorm_Unrelated, 'omitnan')/sqrt(length(Qnorm_Unrelated))];
er = errorbar(X,Y,Y_se/2,Y_se/2);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('TV-Modularity')
ylabel('Euclidean distance')
axis square
hold off

% export_fig(figure(1), fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus_barplot.tif']), '-r300')
exportgraphics(figure(1), fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus_barplot.tif']), 'Resolution', 300)

%% Visual Inspection of the variables (DMNCON, DMNFPN)

figure(2)
clf
set(figure(2), 'color', 'white', 'defaultAxesFontSize', 20)
tiledlayout(1, 2)

nexttile
Y = [mean(DMNDANnorm_MZ, 'omitnan'), mean(DMNDANnorm_DZ, 'omitnan'), mean(DMNDANnorm_NT, 'omitnan'), mean(DMNDANnorm_Unrelated, 'omitnan')];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [std(DMNDANnorm_MZ, 'omitnan')/sqrt(length(DMNDANnorm_MZ)), std(DMNDANnorm_DZ, 'omitnan')/sqrt(length(DMNDANnorm_DZ)), std(DMNDANnorm_NT, 'omitnan')/sqrt(length(DMNDANnorm_NT)), std(DMNDANnorm_Unrelated, 'omitnan')/sqrt(length(DMNDANnorm_Unrelated))];
er = errorbar(X,Y,Y_se/2,Y_se/2);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('TV-FC DMN-DAN')
ylabel('Euclidean distance')
axis square
hold off

nexttile
Y = [mean(DMNCONnorm_MZ, 'omitnan'), mean(DMNCONnorm_DZ, 'omitnan'), mean(DMNCONnorm_NT, 'omitnan'), mean(DMNCONnorm_Unrelated, 'omitnan')];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [std(DMNCONnorm_MZ, 'omitnan')/sqrt(length(DMNCONnorm_MZ)), std(DMNCONnorm_DZ, 'omitnan')/sqrt(length(DMNCONnorm_DZ)), std(DMNCONnorm_NT, 'omitnan')/sqrt(length(DMNCONnorm_NT)), std(DMNCONnorm_Unrelated, 'omitnan')/sqrt(length(DMNCONnorm_Unrelated))];
er = errorbar(X,Y,Y_se/2,Y_se/2);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('TV-FC DMN-CON')
ylabel('Euclidean distance')
axis square
hold off

exportgraphics(figure(2), fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus_barplot_exploratory.tif']), 'Resolution', 300)


%% Run one-way ANOVA of sibling status
delete(findall(0))
for m = 1:length(testing_metrics)
    [anova_p.(testing_metrics{m}),anova_tab.(testing_metrics{m}) ,anova_stats.(testing_metrics{m})] = anova1(eval(['sibling_pair_' testing_metrics{m}]), sibling_pair_label, 'off');
    title(['One-way ANOVA of sibling status: ' testing_metrics{m}])
    close(figure(1))
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus.mat']), 'anova_*', 'sibling_pair*', 'TPnorm*', 'FOnorm*', 'Qnorm*', 'FCnorm*')

%% Run one-way ANOVA of sibling status with age covariate
% delete(findall(0))
age_diff = vertcat(age_MZ(:, 1)-age_MZ(:, 2), age_DZ(:, 1)-age_DZ(:, 2), age_NT(:, 1)-age_NT(:, 2), age_Unrelated(:, 1)-age_Unrelated(:, 2));
age_diff = abs(age_diff);

for m = 1:length(testing_metrics)
    [anovan_p.(testing_metrics{m}),anovan_tab.(testing_metrics{m}) ,anovan_stats.(testing_metrics{m})] = anovan(eval(['sibling_pair_' testing_metrics{m}]), {sibling_pair_label, age_diff}, 'Continuous', 2, 'varnames', {'group', 'age_diff'});
    title({['One-way ANOVA of sibling status: ' testing_metrics{m}]; 'all same sex, age-diff cov'})
    close(figure(1))
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus.mat']), 'anovan_*', 'age_diff', '-append')

%% Run one-way ANOVA of sibling status with age and FD covariates
delete(findall(0))
load(fullfile(dir_home, 'scripts2_hmm_main', 'twinid_FD.mat'), ['FD_*_K' num2str(K)]) %not resampled
FD_MZ = eval(['FD_MZ_K' num2str(K)]);
FD_DZ = eval(['FD_DZ_K' num2str(K)]);
FD_NT = eval(['FD_NT_K' num2str(K)]);
FD_Unrelated = eval(['FD_Unrelated_K' num2str(K)]);

FD_diff = vertcat(FD_MZ(:, 1)-FD_MZ(:, 2), FD_DZ(:, 1)-FD_DZ(:, 2), FD_NT(:, 1)-FD_NT(:, 2), FD_Unrelated(:, 1)-FD_Unrelated(:, 2));
FD_diff = abs(FD_diff);

for m = 1:length(testing_metrics)
    [anovan2_p.(testing_metrics{m}),anovan2_tab.(testing_metrics{m}) ,anovan2_stats.(testing_metrics{m})] = anovan(eval(['sibling_pair_' testing_metrics{m}]), {sibling_pair_label, age_diff, FD_diff}, 'Continuous', [2,3], 'varnames', {'group', 'age_diff', 'FD_diff'});
    title({['One-way ANOVA of sibling status: ' testing_metrics{m}]; 'all same sex, age-diff, FD-diff cov'})
    close(figure(1))
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus.mat']), 'anovan2_*', 'FD_diff', '-append')

%% Report
for m = 1:length(testing_metrics)
    report{m, 1} = testing_metrics{m};
    report{m, 2} = anovan2_tab.(testing_metrics{m}){2, 6};
    report{m, 3} = anovan2_tab.(testing_metrics{m}){2, 7};
    report{m, 4} = anovan2_tab.(testing_metrics{m}){2, 2}/(anovan2_tab.(testing_metrics{m}){2, 2}+anovan2_tab.(testing_metrics{m}){5, 2});
end

%% Run one-way ANOVA of sibling status with age and alternative FD covariates
% delete(findall(0))
clear FD_*
load(fullfile(dir_home, 'scripts2_hmm_main', 'twinid_FD.mat'), 'FD_MZ', 'FD_DZ', 'FD_NT', 'FD_Unrelated') %not resampled

FD_diff2 = vertcat(FD_MZ(:, 1)-FD_MZ(:, 2), FD_DZ(:, 1)-FD_DZ(:, 2), FD_NT(:, 1)-FD_NT(:, 2), FD_Unrelated(:, 1)-FD_Unrelated(:, 2));
FD_diff2 = abs(FD_diff2);

for m = 1:length(testing_metrics)
    [anovan3_p.(testing_metrics{m}),anovan3_tab.(testing_metrics{m}) ,anovan3_stats.(testing_metrics{m})] = anovan(eval(['sibling_pair_' testing_metrics{m}]), {sibling_pair_label, age_diff, FD_diff2}, 'Continuous', [2,3], 'varnames', {'group', 'age_diff', 'altFD_diff'});
    title({['One-way ANOVA of sibling status: ' testing_metrics{m}]; 'all same sex, age-diff, altFD-diff cov'})
    close(figure(1))
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anova1_of_siblingstatus.mat']), 'anovan3_*', 'FD_diff2', '-append')

%% Alternative Report
for m = 1:length(testing_metrics)
    report_alt{m, 1} = testing_metrics{m};
    report_alt{m, 2} = anovan3_tab.(testing_metrics{m}){2, 6};
    report_alt{m, 3} = anovan3_tab.(testing_metrics{m}){2, 7};
    report_alt{m, 4} = anovan3_tab.(testing_metrics{m}){2, 2}/(anovan3_tab.(testing_metrics{m}){2, 2}+anovan3_tab.(testing_metrics{m}){5, 2});
end
