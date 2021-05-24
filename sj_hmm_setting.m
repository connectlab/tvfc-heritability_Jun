%% SETUP THE MATLAB PATHS AND FILE NAMES
% Code used in Vidaurre et al. (2017) PNAS
%
% Detailed documentation and further examples can be found in:
% https://github.com/OHBA-analysis/HMM-MAR
% This pipeline must be adapted to your particular configuration of files.
clear; clc;

dir_scripts = pwd; 
dir_rawHCP = fullfile(dir_scripts, 'RAW_DATA','HCP_YoungAdult_1200');
answer=inputdlg({'K4 or K6'},'Input condition',1,{'4 or 6'});

dir_node = fullfile(dir_rawHCP, 'HCP1200_Parcellation_Timeseries_Netmats', ['NodeTimeseries_ICAd300_ts2']);
dir_netmat = fullfile(dir_rawHCP, 'HCP1200_Parcellation_Timeseries_Netmats', ['netmats_ICAd300_ts2']);
dir_groupica = fullfile(dir_rawHCP, 'HCP1200_Parcellation_Timeseries_Netmats', 'groupICA_3T_HCP1200_MSMAll', 'groupICA', ['groupICA_3T_HCP1200_MSMAll_d300.ica']);
dir_output = fullfile(dir_scripts, ['HMM_d300_xSubcorticals']); mkdir(dir_output);

%% xSubcorticals
yeo_label = {'Visual', 'SMN', 'DAN', 'CON', 'Limbic', 'FPN', 'DMN', 'Subcorticals'};
load(fullfile(dir_scripts, 'ica_wholelabel_info.mat'), 'parcel_info')
parcel_info = parcel_info(2:end, :);

dICA = 1:300;
ic_selected = strcmpi(parcel_info(:, 9), 'Subcorticals');
ic_selected = dICA(~ic_selected);

if ~exist(fullfile(dir_node, 'xSubcorticals'), 'dir') == 1
    mkdir(fullfile(dir_node, 'xSubcorticals'))
    
    txtfiles = cellstr(ls(fullfile(dir_node, '*.txt')));
    NSub = length(txtfiles); % no. subjects
    
    for j=1:NSub
        temp_load = load(fullfile(dir_node, txtfiles{j}));
        temp_load = temp_load(:, ic_selected);
        eval(['X' txtfiles{j}(1:end-4) '=temp_load;']);
        save(fullfile(dir_node, 'xSubcorticals', [txtfiles{j}(1:end-4) '.mat']), ['X' txtfiles{j}(1:end-4)])
        clear temp_load
        clear(['X' txtfiles{j}(1:end-4)])
    end
    clear txtfiles NSub j;
end

matfiles = cellstr(ls(fullfile(dir_node, 'xSubcorticals', '*.mat')));
NSub = length(matfiles); % no. subjects
TimePoints = 4800;

f = cell(NSub,1);
T = cell(NSub,1);

for j=1:NSub
    f{j} = fullfile(dir_node, 'xSubcorticals', matfiles{j});
    % load(f{j});
    T{j} = [1200 1200 1200 1200];
end

%----------- selecting ICs & Yeo definition
parcel_info = parcel_info(ic_selected, 9);

dICA_new = 1:length(ic_selected);
DMN = dICA_new(contains(parcel_info, 'DMN'));
FPN = dICA_new(contains(parcel_info, 'FPN'));
DAN = dICA_new(contains(parcel_info, 'DAN'));
CON = dICA_new(contains(parcel_info, 'CON'));
SMN = dICA_new(contains(parcel_info, 'SMN'));
Visual = dICA_new(contains(parcel_info, 'Visual'));
Limbic = dICA_new(contains(parcel_info, 'Limbic'));

IC_reorder = [DMN, FPN, DAN, CON, SMN, Visual, Limbic];
IC_label = horzcat(repmat(yeo_label(7), 1, length(DMN)),...
    repmat(yeo_label(6), 1, length(FPN)),...
    repmat(yeo_label(3), 1, length(DAN)),...
    repmat(yeo_label(4), 1, length(CON)),...
    repmat(yeo_label(2), 1, length(SMN)),...
    repmat(yeo_label(1), 1, length(Visual)),...
    repmat(yeo_label(5), 1, length(Limbic)));
IC_IND = horzcat(repmat(1, 1, length(DMN)),...
    repmat(2, 1, length(FPN)),...
    repmat(3, 1, length(DAN)),... 
    repmat(4, 1, length(CON)),... 
    repmat(5, 1, length(SMN)),... 
    repmat(6, 1, length(Visual)),... 
    repmat(7, 1, length(Limbic)));




