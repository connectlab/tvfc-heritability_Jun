run sj_hmm_setting

load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

%% Entire subject pool-based metrics: FC matrix (with selective ICs), TP, FO.
%% 1. Calculating subject-wise FC matrix without normalization (nets_netmats)

ts = nets_load(fullfile(dir_node, 'xSubcorticals'), 0.72,1,1);
netmats1=nets_netmats(ts,0,'corr');

for s = 1:NSub
    netmats1_sub{s} = reshape(netmats1(s, :), length(IC_reorder), length(IC_reorder));
    netmats1_sub_fischer_z{s} = atanh(netmats1_sub{s}(IC_reorder, IC_reorder)); %Fisher r-z
    netmats1_z(:, :, s) = netmats1_sub_fischer_z{s};
end

Mnet1_xSubcorticals = mean(netmats1_z, 3);
% [Znet_xSubcorticals,Mnet1_xSubcorticals]=nets_groupmean(netmats1_z,1,1);

save(fullfile(dir_output, 'netmats1.mat'), 'ts', 'netmats1*','Mnet*');%'Znet*', 

%% 2. Transition Proability
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '.mat']))

figure(r)
clf
set(figure(r),'color','white') % Make the background white:

GammaSessMean = squeeze(mean(reshape(Gamma,[1200 4 NSub K]),1));
GammaSubMean = squeeze(mean(GammaSessMean,1));
[~,pca1] = pca(GammaSubMean','NumComponents',1);
[~,ord] = sort(pca1);

% Figure 2A Transition probability Matrix
figs(1) = subplot(1,2,1);
P = getTransProbs(hmm);
% P = hmm.P;
%P(eye(K)==1) = 0;
% for j=1:K, P(j,:) = P(j,:) / sum(P(j,:)); end
figs_matrix{1} = P(ord,ord);
figs_matrix{1}(1:(size(figs_matrix{1},1)+1):end) = NaN;
% by SJ
imAlpha{1} = ones(size(figs_matrix{1}));
imAlpha{1}(isnan(figs_matrix{1}))=0;
imagesc(figs_matrix{1},'AlphaData',imAlpha{1}, [0 round(max(max(P)),1)]);
set(figs(1),'color','white') % Make the background white:
colorbar
colormap('jet')
axis square
hold on
for j=0:K+1
    plot([0 K+1] - 0.5,[j j] + 0.5,'k','LineWidth',1)
    plot([j j] + 0.5,[0 K+1] - 0.5,'k','LineWidth',1)
end
title(figs(1),'Transition probability Matrix','fontsize',10);


% Figure 2B FO probability matrix
figs(2) = subplot(1,2,2);
figs_matrix{2} = corr(GammaSubMean(:,ord));
figs_matrix{2}(1:(size(figs_matrix{2},1)+1):end) = NaN;
% by SJ
imAlpha{2} = ones(size(figs_matrix{2}));
imAlpha{2}(isnan(figs_matrix{2}))=0;
imagesc(figs_matrix{2},'AlphaData',imAlpha{2}, [-1.0 1.0]);
colormap(figs(2),jet)
set(figs(2),'color','white') % Make the background white:
colorbar
colormap('jet')
axis square
hold on
for j=0:K+1
    plot([0 K+1] - 0.5,[j j] + 0.5,'k','LineWidth',1)
    plot([j j] + 0.5,[0 K+1] - 0.5,'k','LineWidth',1)
end
title(figs(2),'Fractional Occupancy matrix','fontsize',10);

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_ord.mat']), 'ord');
export_fig(figure(r), fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_fig_TPFO.tif']), '-r300');

TP_ord = P(ord,ord);

% vectorize TP
ct = 0;
for i = 1:K
    for j = 1:K
        if i == j, else ct = ct+1; TP_ord_vector(ct) = TP_ord(i, j); TP_coord(ct, :) = [i, j];  end
    end
end

TP_ord_vector_z = normalize(TP_ord_vector, 'zscore'); % normalize on vectorized TP

% reshape to matrix
TP_ord_z = zeros(K);
for i = 1:length(TP_coord)
    TP_ord_z(TP_coord(i, 1), TP_coord(i, 2)) = TP_ord_vector_z(i);
end

%% 3. Fractional Occupancy
FO = corr(GammaSubMean(:,ord));

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_entiresubjs.mat']),'TP_ord*', 'TP_coord', 'FO');



