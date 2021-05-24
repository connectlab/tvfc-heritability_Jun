run sj_hmm_setting

%% Set up the HMM parameters
K = str2double(answer{2}); % no. states Increasing K above 6 did not change the topologies of the most prominent states or their task profiles. (Quinn et al., 2018)
r_start = 1;% Quinn (2018) Increasing K above 6 did not change the topologies of the most prominent states or their task profiles.

repetitions = 5; % to run it multiple times (keeping all the results)
TR = 0.72;
use_stochastic = 1; % set to 1 if you have loads of data

options = struct();
options.K = K; % number of states
options.order = 0; % no autoregressive components
options.zeromean = 0; % model the mean
options.covtype = 'full'; % full covariance matrix
options.Fs = 1/TR;
options.verbose = 1;
options.standardise = 1;
options.inittype = 'HMM-MAR';
options.cyc = 500;
options.initcyc = 10;
options.initrep = 3;

% stochastic options
if use_stochastic
    options.BIGNbatch = round(NSub/30);
    options.BIGtol = 1e-7;
    options.BIGcyc = 500;
    options.BIGundertol_tostop = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

%% We run the HMM multiple times
for r = 1:repetitions
    [hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist, feterms, rho]= hmmmar(f,T,options);
    save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '.mat']),'hmm', 'Gamma', 'Xi', 'vpath', 'GammaInit', 'residuals', 'fehist', 'feterms', 'rho')
    disp(['RUN ' num2str(r)])
end

for r = 1:repetitions
    load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '.mat']),'hmm')
    
    [Gamma,Xi] = hmmdecode(f,T,hmm,0);
    
    for i = 1:K        
        eval(['Xi_dim3_' num2str(i) '=Xi(:, :, i);'])
        if exist(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']), 'file') 
            save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']),['Xi_dim3_' num2str(i)], '-append')
        else
           save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']),['Xi_dim3_' num2str(i)]) 
        end        
        clear(['Xi_dim3_' num2str(i)])        
    end
    
    clear Gamma Xi hmm
end







