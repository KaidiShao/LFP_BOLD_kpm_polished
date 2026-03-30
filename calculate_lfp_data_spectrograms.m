% simulation of applying various KMD methods (mainly EDMD) to Nikos's LFP signals
% EDMD: https://venturi.soe.ucsc.edu/sites/default/files/Koopman_0.pdf
% author: Kaidi Shao (20230708)

clear;
restoredefaultpath;

onedrive_path = 'D:\Onedrive\'; %%%%%%%% for ICPBR workstation %%%%%%%%%
% onedrive_path = 'D:\Kaidi\Onedrive\';%%%%%%%% MPI KYB DESKTOP %%%%%%%%
% onedrive_path = 'C:\Users\skd\OneDrive\'; %%%%%%%% DELL LAPTOP %%%%%%%%

% project/code paths
proj_path_main = strcat(onedrive_path,'\ICPBR\Alberta\koopman_events\pipeline_v2\');
addpath(genpath(proj_path_main));
% proj_path_sub = strcat(onedrive_path,'\updated-desnap-with-causality\');
% addpath(genpath(proj_path_sub));

% data paths
data_input_path = 'D:\DataPons\';
data_output_path = 'D:\DataPons_processed\';

% utility function's path
addpath(genpath(strcat(onedrive_path,'\util_functions\')));
addpath(genpath(strcat(onedrive_path,'\Toolbox\eeglab10_2_5_8b\')));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(strcat(data_input_path, ExpSettings.DataFileName, '.mat'));

regions = ExpSettings.RegionAll; % procedure before EDMD consider all channels

lfp_mat = [];
N_dim_region = zeros(1, size(regions, 2));
idx_region = [];
for n_region = 1:size(regions, 2)
    N_dim_region(n_region) = size(eval(strcat(regions{n_region},'_lfp'))', 1);
    lfp_mat = [lfp_mat; eval(strcat(regions{n_region},'_lfp'))'];
    idx_region = [idx_region n_region*ones(1,N_dim_region(n_region))];
end
% X_sr = lfp_mat(1:end-1,:)';
% Y_sr = lfp_mat(2:end,:)';
state_var = lfp_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check saving folders %%%%%%%%%%%%%%%%%%%%%%%%
FigSavePath = strcat(ExpSettings.ExpPath, '\figures\', ExpSettings.DataFileName,'\',ExpSettings.DataType,'\',ExpSettings.KpmMethod,'_',ExpSettings.DictType);
if exist(FigSavePath, 'dir')~=7
    mkdir(FigSavePath)
end
    cd(FigSavePath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% perform classical event detection precedure %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params.Fs = 667;
% Params.passband1 = [5,15]; % according to Nikos's Nature paper
% 
% B1 = fir1(49, Params.passband1/(.5*Params.Fs), 'bandpass');
% D = filter(B1,1,state_var(:,:)')'; % forward filtering
% D(:, 1:24) = []; % remove FIR-induced delays 
% D(:, size(D,2)+[1:24]) = nan; % remove FIR-induced delays 
% 
% d0 = nan(1, size(state_var, 1));
% loc = cell(1, size(state_var, 1));
% loc_peak = cell(1, size(state_var, 1));
% 
% for n_dim = 1:size(state_var, 1)
%     d0(n_dim) = nanmean(D(n_dim,:))+3.5*nanstd(D(n_dim, :));
%     loc{n_dim} = find(D(n_dim, :)>d0(n_dim)); loc{n_dim}(1)=[];
%     loc_t = D(n_dim, :)>d0(n_dim);
% 
%     L_start = 50; L_extract = 100;
%     loc_peak{n_dim} = find_peak_loc(D(n_dim,:), loc{n_dim}, L_extract);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% calculate spectrogram observables %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N_dim, N_sample] = size(state_var);
N_region = numel(regions);
Fs = 667;

Nfreqs = 251;
N_seg = 10;
N_sample_seg = floor(N_sample / N_seg);

spectrogram_save_path = fullfile(data_output_path, 'spectrograms');
if exist(spectrogram_save_path, 'dir') ~= 7
    mkdir(spectrogram_save_path);
end

file_mean_abs = fullfile(spectrogram_save_path, ...
    [ExpSettings.DataName, '_mean_spectrograms_abs.mat']);
file_mean_complex = fullfile(spectrogram_save_path, ...
    [ExpSettings.DataName, '_mean_spectrograms_complex.mat']);

tmpall_mean_abs = nan(Nfreqs, N_sample, N_region);
tmpall_mean_complex = complex( ...
    nan(Nfreqs, N_sample, N_region), ...
    nan(Nfreqs, N_sample, N_region));

if exist(file_mean_abs, 'file') == 2 && exist(file_mean_complex, 'file') == 2
    S_abs = load(file_mean_abs, 'tmpall_mean_abs', 'freqs', 'timesout', 'regions');
    S_complex = load(file_mean_complex, 'tmpall_mean_complex');

    tmpall_mean_abs = S_abs.tmpall_mean_abs;
    tmpall_mean_complex = S_complex.tmpall_mean_complex;
    freqs = S_abs.freqs;
    timesout = S_abs.timesout;
    regions = S_abs.regions;
else
    for n_seg = 1:N_seg
        seg_start = N_sample_seg * (n_seg - 1) + 1;

        if n_seg < N_seg
            seg_end = N_sample_seg * n_seg;
        else
            seg_end = N_sample;
        end

        t_seg = seg_start:seg_end;

        tmpall = complex( ...
            nan(Nfreqs, numel(t_seg), N_dim), ...
            nan(Nfreqs, numel(t_seg), N_dim));

        for n_dim = 1:N_dim
            [tmpall(:,:,n_dim), freqs, timesout] = timefreqMB( ...
                state_var(n_dim, t_seg)', Fs, 'freqs', [0.1, 250], 'nfreqs', Nfreqs);
        end

        file_seg = fullfile(spectrogram_save_path, ...
            [ExpSettings.DataName, '_spectrograms_seg_', int2str(n_seg), '.mat']);
        save(file_seg, 'tmpall', 'freqs', 'timesout', 'regions', '-v7.3');

        for n_region = 1:N_region
            region_mask = (idx_region == n_region);

            tmpall_mean_complex(:, t_seg, n_region) = ...
                squeeze(mean(tmpall(:,:,region_mask), 3));

            tmpall_mean_abs(:, t_seg, n_region) = ...
                squeeze(mean(abs(tmpall(:,:,region_mask)), 3));
        end
    end

    save(file_mean_abs, 'tmpall_mean_abs', 'freqs', 'timesout', 'regions', '-v7.3');
    save(file_mean_complex, 'tmpall_mean_complex', 'freqs', 'timesout', 'regions', '-v7.3');
end