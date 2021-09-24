run sj_hmm_setting
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

load twinid

subid_MZ = pair_MZ(:);
subid_DZ = pair_DZ(:);
subid_NT = pair_NT(:);
subid_Unrelated = pair_Unrelated(:);

%% Prep for Two-way ANOVA analysis of sibling similarity x brain states.
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'Yeo_NewmanQ', 'DMNFPN')
load(fullfile(dir_hmm, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'FO_SubState', 'TP_Sub_ord_vector')

testing_metrics = {'Q','FC','FO','TP'};


% Redefine the tangent point: subtract surrogate mean from all subjects.
dir_simhmmmar_stats = replace(dir_scripts, 'scripts2_hmm_main_stats', 'scripts5_hmm_simhmmmar_stats');
load(fullfile(dir_simhmmmar_stats, ['HMMrun_K' num2str(K) '_surrogate_mean.mat']))

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
load(fullfile(dir_home, 'scripts2_hmm_main', 'twinid_FD.mat'), ['FD_*_K' num2str(K)]) %not resampled
FD_MZ = eval(['FD_MZ_K' num2str(K)]);
FD_DZ = eval(['FD_DZ_K' num2str(K)]);
FD_NT = eval(['FD_NT_K' num2str(K)]);
FD_Unrelated = eval(['FD_Unrelated_K' num2str(K)]);

FD_diff = vertcat(FD_MZ(:, 1)-FD_MZ(:, 2), FD_DZ(:, 1)-FD_DZ(:, 2), FD_NT(:, 1)-FD_NT(:, 2), FD_Unrelated(:, 1)-FD_Unrelated(:, 2));
FD_diff = abs(FD_diff);

age_diff = vertcat(age_MZ(:, 1)-age_MZ(:, 2), age_DZ(:, 1)-age_DZ(:, 2), age_NT(:, 1)-age_NT(:, 2), age_Unrelated(:, 1)-age_Unrelated(:, 2));
age_diff = abs(age_diff);

for m = 1:length(testing_metrics)
    FD_cov = repmat(FD_diff, size(eval(testing_metrics{m}), 2), 1);
    age_cov = repmat(age_diff, size(eval(testing_metrics{m}), 2), 1);
    data = eval(['sibling_pair_' testing_metrics{m}]);
    sibling_label = eval(['sibling_label_' testing_metrics{m}]);
    state_label = eval(['state_label_' testing_metrics{m}]);
    
    [anovan_p.(testing_metrics{m}),anovan_tab.(testing_metrics{m}),anovan_stats.(testing_metrics{m})] = anovan(data, {sibling_label, state_label, age_cov, FD_cov},'Continuous', [3,4],'varnames',{'sibling status','State K', 'age_diff', 'FD_diff'});
%      'model','interaction',
    %     multcompare(anovan_stats.(testing_metrics{m})); 
     
    clear data sibling_label state_label age_cov FD_cov
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anovan_siblingsimilarity_x_state.mat']), 'anovan_*', 'sibling_*', '*label*')

delete(findall(0)) %close all open windows

%% Report
for m = 1:length(testing_metrics)
    report{m, 1} = testing_metrics{m};
    report{m, 2} = anovan_tab.(testing_metrics{m}){2, 6};
    report{m, 3} = anovan_tab.(testing_metrics{m}){2, 7};
    report{m, 4} = anovan_tab.(testing_metrics{m}){2, 2}/(anovan_tab.(testing_metrics{m}){2, 2}+anovan_tab.(testing_metrics{m}){6, 2});
end
%% Visual Inspection of the mean of the sibling groups across the variables

figure(1), clf
set(figure(1), 'color', 'white','windowstate', 'maximized', 'defaultAxesFontSize', 22)
tiledlayout(1, 4)

pair_FC = reshape(sibling_pair_FC, length(sibling_pair_FC)/K, K);
pair_Q = reshape(sibling_pair_Q, length(sibling_pair_Q)/K, K);
pair_FO = reshape(sibling_pair_FO, length(sibling_pair_FO)/K, K);
pair_TP = reshape(sibling_pair_TP, length(sibling_pair_TP)/(K*(K-1)), K*(K-1));

%---------------------------------------
nexttile
clear X Y Y_se er 
Y = [mean(mean(pair_FO(1:length(pair_MZ), :), 'omitnan')),...
    mean(mean(pair_FO(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')),...
    mean(mean(pair_FO(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')),...
    mean(mean(pair_FO(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan'))];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [mean(std(pair_FO(1:length(pair_MZ), :), 'omitnan')/sqrt(length(subid_MZ))),...
    mean(std(pair_FO(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')/sqrt(length(subid_DZ))),...
    mean(std(pair_FO(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')/sqrt(length(subid_NT))),...
    mean(std(pair_FO(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan')/sqrt(length(subid_Unrelated)))];
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Fractional Occupancy')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

%---------------------------------------
nexttile
clear X Y Y_se er 
Y = [mean(mean(pair_TP(1:length(pair_MZ), :), 'omitnan')),...
    mean(mean(pair_TP(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')),...
    mean(mean(pair_TP(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')),...
    mean(mean(pair_TP(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan'))];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [mean(std(pair_TP(1:length(pair_MZ), :), 'omitnan')/sqrt(length(subid_MZ))),...
    mean(std(pair_TP(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')/sqrt(length(subid_DZ))),...
    mean(std(pair_TP(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')/sqrt(length(subid_NT))),...
    mean(std(pair_TP(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan')/sqrt(length(subid_Unrelated)))];
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Transition Probability')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

%---------------------------------------
nexttile
clear X Y Y_se er 
Y = [mean(mean(pair_FC(1:length(pair_MZ), :), 'omitnan')),...
    mean(mean(pair_FC(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')),...
    mean(mean(pair_FC(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')),...
    mean(mean(pair_FC(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan'))];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [mean(std(pair_FC(1:length(pair_MZ), :), 'omitnan')/sqrt(length(subid_MZ))),...
    mean(std(pair_FC(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')/sqrt(length(subid_DZ))),...
    mean(std(pair_FC(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')/sqrt(length(subid_NT))),...
    mean(std(pair_FC(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan')/sqrt(length(subid_Unrelated)))];
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('TV-FC DMN-FPN')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

%---------------------------------------
nexttile
clear X Y Y_se er 
Y = [mean(mean(pair_Q(1:length(pair_MZ), :), 'omitnan')),...
    mean(mean(pair_Q(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')),...
    mean(mean(pair_Q(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')),...
    mean(mean(pair_Q(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan'))];
X = categorical({'MZ','DZ','NT','Unrelated'});
X = reordercats(X,{'MZ','DZ','NT','Unrelated'});
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
Y_se = [mean(std(pair_Q(1:length(pair_MZ), :), 'omitnan')/sqrt(length(subid_MZ))),...
    mean(std(pair_Q(1+length(pair_MZ):length(pair_MZ)+length(pair_DZ), :), 'omitnan')/sqrt(length(subid_DZ))),...
    mean(std(pair_Q(1+length(pair_MZ)+length(pair_DZ):length(pair_MZ)+length(pair_DZ)+length(pair_NT), :), 'omitnan')/sqrt(length(subid_NT))),...
    mean(std(pair_Q(1+length(pair_MZ)+length(pair_DZ)+length(pair_NT):length(pair_MZ)+length(pair_DZ)+length(pair_NT)+length(pair_Unrelated), :), 'omitnan')/sqrt(length(subid_Unrelated)))];
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('TV-Modularity')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

exportgraphics(figure(1), fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anovan_siblingsimilarity_x_state_barplot1.tif']), 'Resolution', 300)

close(figure(1))


%% Visual Inspection of the mean of the variables across the sibling groups

figure(2), clf
set(figure(2), 'color', 'white','windowstate', 'maximized', 'defaultAxesFontSize', 22)
tiledlayout(1, 4)

pair_FC = reshape(sibling_pair_FC, length(sibling_pair_FC)/K, K);
pair_Q = reshape(sibling_pair_Q, length(sibling_pair_Q)/K, K);
pair_FO = reshape(sibling_pair_FO, length(sibling_pair_FO)/K, K);
pair_TP = reshape(sibling_pair_TP, length(sibling_pair_TP)/(K*(K-1)), K*(K-1));

%---------------------------------------
nexttile
clear X Y Y_se er Z
for kk = 1:K
    Z{kk} = ['K' num2str(kk)];
end
X = reordercats(categorical(Z),Z);
Y = mean(pair_FO, 'omitnan');
for kk = 1:K
    Y_se(kk) = std(pair_FO(:,kk), 'omitnan')/sqrt(length(pair_FO(~isnan(pair_FO(:, kk)), kk)));
end
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Fractional Occupancy')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

%---------------------------------------
nexttile
clear X Y Y_se er Z
for kk = 1:K*(K-1)
    Z{kk} = ['TPcomp' num2str(kk)];
end
X = reordercats(categorical(Z),Z);
Y = mean(pair_TP, 'omitnan');
for kk = 1:K*(K-1)
    Y_se(kk) = std(pair_TP(:,kk), 'omitnan')/sqrt(length(pair_TP(~isnan(pair_TP(:, kk)), kk)));
end
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
% xtickangle(90)
% xlabel(X, 'FontSize', 12)
hold on
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('Transition Probability')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

%---------------------------------------
nexttile
clear X Y Y_se er Z
for kk = 1:K
    Z{kk} = ['K' num2str(kk)];
end
X = reordercats(categorical(Z),Z);
Y = mean(pair_FC, 'omitnan');
for kk = 1:K
    Y_se(kk) = std(pair_FC(:,kk), 'omitnan')/sqrt(length(pair_FC(~isnan(pair_FC(:, kk)), kk)));
end
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('TV-FC DMN-FPN')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er

%---------------------------------------
nexttile
clear X Y Y_se er Z
for kk = 1:K
    Z{kk} = ['K' num2str(kk)];
end
X = reordercats(categorical(Z),Z);
Y = mean(pair_Q, 'omitnan');
for kk = 1:K
    Y_se(kk) = std(pair_Q(:,kk), 'omitnan')/sqrt(length(pair_Q(~isnan(pair_Q(:, kk)), kk)));
end
bar(X, Y, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1)
hold on
er = errorbar(X,Y,Y_se/2,Y_se/2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title('TV-Modularity')
ylabel('Euclidean distance')
axis square
hold off
clear X Y er


exportgraphics(figure(2), fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_anovan_siblingsimilarity_x_state_barplot2.tif']), 'Resolution', 300)

delete(findall(0))