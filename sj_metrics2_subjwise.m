run sj_hmm_setting

load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_free_energy.mat']), 'main_inference_run')
r = main_inference_run;

%% Subject-wise measures (independent of state)
%% Newman's Modularity of each subject with predefined Ci
load(fullfile(dir_output, 'netmats1.mat'), 'netmats1_sub_fischer_z');

for s = 1:NSub
    Mnet1_xSubcorticals{s} = netmats1_sub_fischer_z{s}; %Reordered already.
    
    gamma = 1;                                    %resolution parameter
    Ci=IC_IND';                                   %pre-defined community indices_SJ
    n_vertices=length(Mnet1_xSubcorticals{s});       %number of vertices
    d=sum(Mnet1_xSubcorticals{s});                   %degree
    n_edges=sum(d);                               %number of edges (each undirected edge is counted twice)
    B=Mnet1_xSubcorticals{s}-gamma*(d.'*d)/n_edges;  %modularity matrix
    m=Ci(:,ones(1,n_vertices));                   %compute modularity
    Q=~(m-m.').*B/n_edges;
    Q=sum(Q(:));
    
    Yeo_NewmanQ(s, 1) = Q;
    
    % DMN x TPN mean functional connectivity value
    DMNFPN(s, 1) = mean(mean(Mnet1_xSubcorticals{s}(IC_IND==1, IC_IND==2)));
    DMNDAN(s, 1) = mean(mean(Mnet1_xSubcorticals{s}(IC_IND==1, IC_IND==3)));
    DMNCON(s, 1) = mean(mean(Mnet1_xSubcorticals{s}(IC_IND==1, IC_IND==4)));
    
    clear gamma Ci n_vertices d n_edges B m Q
end
save(fullfile(dir_output, 'metrics_subjwise_corticalFCM.mat'),'Mnet*', 'Yeo_NewmanQ', 'DMM*N');


%% State-wise measures
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '.mat']));

for s = 1:length(matfiles)
    load(fullfile(dir_node, 'xSubcorticals', matfiles{s}));
    subj = eval(['X' matfiles{s}(1:6)]);
    subj_vpath = vpath(TimePoints*(s-1)+1:TimePoints*s, 1);
    
    %% reconstruct the state-wise functional connectivity matrix of a given subject
    for i = 1:K
        subj_state = subj(subj_vpath==i, :);
        
        if size(subj_state, 1) < 2
            Yeo_NewmanQ(s, i) = NaN;
            DMNFPN(s, i) = NaN;
            DMNDAN(s, i) = NaN;
            DMNCON(s, i) = NaN;
            
        else
            subj_Mnet1 = corr(subj_state);
            subj_Mnet1 = subj_Mnet1(IC_reorder, IC_reorder);
            subj_Mnet1(eye(size(subj_Mnet1, 1))>0)=0;
            subj_Mnet1_z = atanh(subj_Mnet1); % reordered
            %             if sum(isinf(subj_Mnet1_z), 'all')
            
            Mnet1_xSubcorticals{s, i} = subj_Mnet1_z; %To be consistent with sj_metrics1_entiresubjs.m
            
            %% Newman's Modularity of the entire subjects with predefined Ci
            % Following code was modified from modularity_und.m of the BCT.
            % BCT(https://sites.google.com/site/bctnet/measures/list#TOC-Clustering-and-Community-Structure)
            
            gamma = 1;                                    %resolution parameter
            Ci=IC_IND';                                   %pre-defined community indices_SJ
            n_vertices=length(Mnet1_xSubcorticals{s, i});       %number of vertices
            d=sum(Mnet1_xSubcorticals{s, i});                   %degree
            n_edges=sum(d);                               %number of edges (each undirected edge is counted twice)
            B=Mnet1_xSubcorticals{s, i}-gamma*(d.'*d)/n_edges;  %modularity matrix
            m=Ci(:,ones(1,n_vertices));                   %compute modularity
            Q=~(m-m.').*B/n_edges;
            Q=sum(Q(:));
            
            if Q <0
                Yeo_NewmanQ(s, i) = NaN;
            else
                Yeo_NewmanQ(s, i) = Q;
            end
            
            clear gamma Ci n_vertices d n_edges B m Q
            
            %% compute DMN x TPN mean functional connectivity
            DMNFPN(s, i) = mean(mean(Mnet1_xSubcorticals{s, i}(IC_IND==1, IC_IND==2)));
            DMNDAN(s, i) = mean(mean(Mnet1_xSubcorticals{s, i}(IC_IND==1, IC_IND==3)));
            DMNCON(s, i) = mean(mean(Mnet1_xSubcorticals{s, i}(IC_IND==1, IC_IND==4)));
            
        end
        clear subj_state subj_Mnet1*
    end
    clear ('subj*', ['X' matfiles{s}(1:6)])
end
save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise_corticalFCM.mat']), 'Yeo_NewmanQ', 'DMN*', 'Mnet1_xSubcorticals');%'DMNFPN'
clear YeoQ_state_timeseries DMN*

%% Fractional Occupancy & Switching Rate
% Fractional Occupancy for each state for each subject
% if dim = 2, how much time each subject/trial/session spends in each state (i.e. the average state probability across time, per session or subject)
% Switching Rate: a measure of the state switching rate for each session/subject, and can be understood as a measure of stability per subject.

dim = 2;
[FO_SubState, ntrials] = getFractionalOccupancy(Gamma,T,options,dim);

switchingRate = getSwitchingRate(Gamma,T,options);

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'FO_SubState*', 'switchingRate'); %'FO_max',
clear dim FO_* switchingRate

%% Transition Probability for each subject_SJ
% Following was modifed from getMaskedTransProbMats.m

% Masks = {1:TimePoints};
for s = 1:NSub
    Masks{s} = [TimePoints*(s-1)+1:TimePoints*s];
    Masks_Xi{s} = [(TimePoints-4)*(s-1)+1:(TimePoints-4)*s];
    T2(s, :) = [TimePoints*(s-1)+1200, 1200, 1200, 1200];
end

load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '.mat']), 'hmm', 'Gamma');
load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_ord.mat']));

for i = 1:K % [Gamma2,Xi] = hmmdecode(f,T,hmm,0);
    load(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_Xi.mat']), ['Xi_dim3_' num2str(i)]);
    Xi(:, :, i) = eval(['Xi_dim3_' num2str(i)]);
    clear(['Xi_dim3_' num2str(i)])
end

order = hmm.train.maxorder;
embeddedlags = abs(hmm.train.embeddedlags);
L = order + embeddedlags(1) + embeddedlags(end);

if ~iscell(Masks), Masks = {Masks}; end

np = length(Masks);
% P = cell(1,np); Pi = cell(1,np);
% we do not care about the grouping imposed in the inference
if isfield(hmm.train,'grouping'), hmm.train = rmfield(hmm.train,'grouping'); end
if ~isfield(hmm.train,'Pstructure'), hmm.train.Pstructure = true(hmm.K); end
if ~isfield(hmm.train,'Pistructure'), hmm.train.Pistructure = true(1,hmm.K); end


for im = 1:np
    %fprintf('Mask %d \n',im)
    T2 = [TimePoints*(im-1)+1200, 1200, 1200, 1200];
    T2_xi = [(TimePoints-4)*(im-1)+1199, 1199, 1199, 1199];
    N = length(T2);
    mask = Masks{im};
    mask_xi = Masks_Xi{im};
    T0 = []; Gamma0 = []; Xi0 = [];
    for n = 1:N
        t0 = sum(T2(1:n-1)); t1 = sum(T2(1:n));
        t0_xi = sum(T2_xi(1:n-1)); t1_xi = sum(T2_xi(1:n));
        ind_ix = mask(mask>=t0+1 & mask<=t1); % the ones belonging to this trial
        ind_ixi = mask_xi(mask_xi>=t0_xi+1 & mask_xi<=t1_xi); % the ones belonging to this trial
        if length(ind_ix)<=L, continue; end
        T0 = [T0; length(ind_ix)];
        if order > 0
            ind_ig = ind_ix(ind_ix>=t0+order+1);
            ind_ig = ind_ig - n*order;
        elseif length(embeddedlags) > 1
            ind_ig = ind_ix((ind_ix>=t0+embeddedlags(1)+1) & (ind_ix<=t1-embeddedlags(end))  );
            ind_ig = ind_ig - (n-1)*L - embeddedlags(1);
        else
            ind_ig = ind_ix;
        end
        %             ind_ixi = ind_ig(1:end-1) - (n-1);
        Gamma0 = cat(1,Gamma0,Gamma(ind_ig,:));
        Xi0 = cat(1,Xi0,Xi(ind_ixi,:,:));
    end
    if isempty(Gamma0), error('Invalid mask?'); end
    hmm0 = hsupdate(Xi0,Gamma0,T0,hmm);
    TP_Sub{im} = hmm0.P; TP_Sub_Pi_stochastic(im, :) = hmm0.Pi;
    
    % sorting based on ord
    for j=1:K, TP_Sub{im}(j,j) = 0; TP_Sub{im}(j,:) = TP_Sub{im}(j,:) / sum(TP_Sub{im}(j,:));  end
    TP_Sub_ord{im} = TP_Sub{im}(ord, ord);
end

save(fullfile(dir_output, ['HMMrun_K' num2str(K) '_rep' num2str(r) '_metrics_subj_statewise.mat']), 'TP_Sub_*', '-append');

clear TP_Sub_* hmm Gamma Xi


