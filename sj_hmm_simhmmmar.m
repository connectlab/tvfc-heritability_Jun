run sj_hmm_setting
dir_output = fullfile(dir_scripts, ['HMM_d' answer{1} '_xSubcorticals']); mkdir(dir_output);

%%
K = 4; % number of states
ndim = 139; % number of brain regions
Fs = 1/TR;
T = TimePoints * ones(NSub,1); % number of data points
iter = 50;

for r = 37
    close all
    
%     load(fullfile(dir_mainhmm, 'HMMrun_K4_rep3.mat'), 'hmm'); % our main hmm model.
    load(fullfile(dir_mainhmm, 'HMMrun_K6_rep1.mat'), 'hmm'); % our main hmm model.
    hmmtrue = hmm; clear hmm
    
    hmmtrue.Pi = ones(1,K); %rand(1,K);
    hmmtrue.Pi = hmmtrue.Pi./sum(hmmtrue.Pi);
    
    %%
    hmmtrue.P = rand(K) + 100 * eye(K);
    for j=1:K
        hmmtrue.P(j,:) = hmmtrue.P(j,:) ./ sum(hmmtrue.P(j,:));
    end
    
    [X,T2,Gammatrue] = simhmmmar(T,hmmtrue,[]);
    
    fc_corr = zeros(ndim, ndim);
    for j = 1:NSub
        t = (1:T2(j))+sum(T2(1:j-1));
        fc_corr = fc_corr+corr(X(t, :));
    end
    fc_corr = fc_corr/NSub;
    
    figure(1)
    clf;
    set(figure(1), 'color', 'white')
    imagesc(fc_corr(IC_reorder, IC_reorder)+diag(nan(1,min(size(X)))))
    axis square; colorbar; colormap jet; caxis([-1 1])
    title('FC of simulated data with random TP');
    
    figure(2)
    clf;
    set(figure(2), 'color', 'white')
    imagesc(hmmtrue.P+diag(nan(1,K)))
    axis square; colorbar
    title('hmmtrue.P random structure');
    
    save(fullfile(dir_output, 'simhmmmar_1st_tprand.mat'), 'hmmtrue', 'X')
    export_fig(figure(1), fullfile(dir_output, 'tprand_X_FC.tif'), '-r300')
    export_fig(figure(2), fullfile(dir_output, 'tprand_hmmtrueP.tif'), '-r300')
    %% SAVE DATA
    
    if ~exist(fullfile(dir_scripts, 'xSubcorticals', ['randsim' num2str(r)]))
        mkdir(fullfile(dir_scripts, 'xSubcorticals', ['randsim' num2str(r)]))
    end
    for j=1:NSub
        f{j, r} = fullfile(dir_scripts, 'xSubcorticals', ['randsim' num2str(r)], matfiles{j});
        Tcell{j, 1} = [1200 1200 1200 1200];
    end
    
    for j = 1:NSub
        t = (1:T2(j))+sum(T2(1:j-1));
        eval(['X' matfiles{j}(1:end-4) '=X(t, :);']);
        save(fullfile(dir_scripts, 'xSubcorticals', ['randsim' num2str(r)], matfiles{j}), ['X' matfiles{j}(1:end-4)])
        clear(['X' matfiles{j}(1:end-4)])
    end
    %% RUN REAL HMM
    options.DirichletDiag = 10; % Change this to alter FO.
    
    [hmm, Gamma, Xi, vpath] = hmmmar(X,T,options);
    
    save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '.mat']), 'hmm', 'Gamma', 'vpath')
    
    %% SAVE Xi
    [Gamma,Xi] = hmmdecode(f(:, r),Tcell,hmm,0);
    
    for i = 1:K
        eval(['Xi_dim3_' num2str(i) '=Xi(:, :, i);'])
        if exist(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']), 'file')
            save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']),['Xi_dim3_' num2str(i)], '-append')
        else
            save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']),['Xi_dim3_' num2str(i)])
        end
        clear(['Xi_dim3_' num2str(i)])
    end
    
    clear Gamma Xi hmm Gamma vpath hmmtrue 
    
end