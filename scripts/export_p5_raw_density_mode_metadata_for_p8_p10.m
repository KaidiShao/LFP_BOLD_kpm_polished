function export_p5_raw_density_mode_metadata_for_p8_p10()
%EXPORT_P5_RAW_DENSITY_MODE_METADATA_FOR_P8_P10 Export raw efun metadata for P8/P10 plots.
%
% The P8/P10 peak CSV files contain density_file and density_index, but they do
% not currently include raw efun timescale metadata.  This helper scans current
% P5 raw thresholded density MAT files and exports their mode_metadata tables so
% the Python v2 plotting script can join timescale/frequency information without
% recomputing xcorr.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'functions')));
processed_root = io_project.get_project_processed_root();
datasets = {'e10aw1', 'e10bv1', 'e10fV1', 'e10gb1', 'e10gh1', ...
    'e10gw1', 'f12m01', 'k13m17', 'k13m23'};
output_dir = fullfile(repo_root, 'results', ...
    'p8_p10_parameter_selection_current_rmsenv_adaptive_v2');
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end
output_file = fullfile(output_dir, 'p5_raw_density_mode_metadata.csv');

T_all = table();
for i_ds = 1:numel(datasets)
    dataset = datasets{i_ds};
    raw_root = fullfile(processed_root, dataset, 'pipeline5_raw_thresholded_density');
    if exist(raw_root, 'dir') ~= 7
        continue;
    end
    files = dir(fullfile(raw_root, '**', '*rmsenv_adaptive*.mat'));
    for i_file = 1:numel(files)
        mat_file = fullfile(files(i_file).folder, files(i_file).name);
        try
            S = load(mat_file, 'mode_metadata', 'D');
        catch ME
            warning('Failed to load mode_metadata from %s: %s', mat_file, ME.message);
            continue;
        end
        if ~isfield(S, 'mode_metadata') || isempty(S.mode_metadata)
            continue;
        end
        M = S.mode_metadata;
        if ~istable(M)
            warning('mode_metadata is not a table in %s.', mat_file);
            continue;
        end
        M = local_add_preferred_timescale_columns(M, S);
        n = height(M);
        prefix = table();
        prefix.dataset = repmat(string(dataset), n, 1);
        prefix.density_file = repmat(string(mat_file), n, 1);
        prefix.condition_folder = repmat(string(local_condition_folder(mat_file)), n, 1);
        prefix.density_name_for_p8_p10 = repmat( ...
            string(local_density_name_from_condition(mat_file)), n, 1);
        M = [prefix, M]; %#ok<AGROW>
        T_all = local_vertcat_table_union(T_all, M); %#ok<AGROW>
    end
end

if isempty(T_all)
    T_all = table();
end
writetable(T_all, output_file);
fprintf('Exported %d raw density metadata rows to:\n  %s\n', ...
    height(T_all), output_file);
end


function T = local_vertcat_table_union(A, B)
if isempty(A)
    T = B;
    return;
end
if isempty(B)
    T = A;
    return;
end

vars_a = string(A.Properties.VariableNames);
vars_b = string(B.Properties.VariableNames);
all_vars = unique([vars_a, vars_b], 'stable');
A = local_add_missing_table_vars(A, all_vars, B);
B = local_add_missing_table_vars(B, all_vars, A);
T = [A(:, cellstr(all_vars)); B(:, cellstr(all_vars))];
end


function T = local_add_missing_table_vars(T, all_vars, template)
vars = string(T.Properties.VariableNames);
for i_var = 1:numel(all_vars)
    name = all_vars(i_var);
    if any(vars == name)
        continue;
    end
    n = height(T);
    if ismember(name, string(template.Properties.VariableNames))
        exemplar = template.(char(name));
        if isstring(exemplar)
            T.(char(name)) = strings(n, 1);
        elseif iscellstr(exemplar) || iscell(exemplar)
            T.(char(name)) = repmat({''}, n, 1);
        elseif isnumeric(exemplar)
            T.(char(name)) = nan(n, 1);
        elseif islogical(exemplar)
            T.(char(name)) = false(n, 1);
        else
            T.(char(name)) = strings(n, 1);
        end
    else
        T.(char(name)) = strings(n, 1);
    end
end
end


function condition = local_condition_folder(mat_file)
parts = split(string(mat_file), filesep);
idx = find(parts == "pipeline5_raw_thresholded_density", 1, 'last');
if isempty(idx) || idx + 1 > numel(parts)
    condition = "";
else
    condition = char(parts(idx + 1));
end
end


function density_name = local_density_name_from_condition(mat_file)
condition = lower(string(local_condition_folder(mat_file)));
if contains(condition, "standardize")
    suffix = "_standardize_rmsenv_adaptive";
else
    suffix = "_rmsenv_adaptive";
end
if contains(condition, "complex") || contains(condition, "csplit")
    density_name = "raw_csplit_q070" + suffix;
elseif contains(condition, "abs")
    density_name = "raw_abs_q070" + suffix;
else
    density_name = "raw_unknown_q070" + suffix;
end
end


function M = local_add_preferred_timescale_columns(M, S)
dt = [];
if isfield(S, 'D') && isstruct(S.D) && isfield(S.D, 'input') && ...
        isstruct(S.D.input) && isfield(S.D.input, 'dt') && ~isempty(S.D.input.dt)
    dt = double(S.D.input.dt);
    dt = dt(1);
end

if ismember('timescale_sec', M.Properties.VariableNames) && ...
        ~ismember('timescale_sec_bilinear', M.Properties.VariableNames)
    M.timescale_sec_bilinear = M.timescale_sec;
end
if ismember('frequency_hz', M.Properties.VariableNames) && ...
        ~ismember('frequency_hz_bilinear', M.Properties.VariableNames)
    M.frequency_hz_bilinear = M.frequency_hz;
end

if ismember('lambda_discrete', M.Properties.VariableNames) && ...
        ~isempty(dt) && isfinite(dt) && dt > 0
    lambda = M.lambda_discrete(:);
    rho = abs(lambda);
    theta = abs(angle(lambda));

    tau = nan(height(M), 1);
    valid_decay = isfinite(rho) & rho > 0 & rho < 1;
    tau(valid_decay) = -dt ./ log(rho(valid_decay));
    tau(isfinite(rho) & rho >= 1) = Inf;

    freq = nan(height(M), 1);
    valid_freq = isfinite(theta);
    freq(valid_freq) = theta(valid_freq) ./ (2 * pi * dt);

    M.timescale_sec_discrete_log = tau;
    M.frequency_hz_discrete_angle = freq;
else
    M.timescale_sec_discrete_log = nan(height(M), 1);
    M.frequency_hz_discrete_angle = nan(height(M), 1);
end

tau_preferred = M.timescale_sec_discrete_log;
tau_source = repmat("discrete_log_abs_evalue", height(M), 1);
missing_tau = isnan(tau_preferred);
if ismember('timescale_sec_bilinear', M.Properties.VariableNames)
    tau_preferred(missing_tau) = M.timescale_sec_bilinear(missing_tau);
    tau_source(missing_tau) = "bilinear_realpart";
end
tau_source(isnan(tau_preferred)) = "missing";
M.timescale_sec_preferred = tau_preferred;
M.timescale_source_preferred = tau_source;

freq_preferred = M.frequency_hz_discrete_angle;
freq_source = repmat("discrete_angle", height(M), 1);
missing_freq = isnan(freq_preferred);
if ismember('frequency_hz_bilinear', M.Properties.VariableNames)
    freq_preferred(missing_freq) = M.frequency_hz_bilinear(missing_freq);
    freq_source(missing_freq) = "bilinear_imag";
end
freq_source(isnan(freq_preferred)) = "missing";
M.frequency_hz_preferred = freq_preferred;
M.frequency_source_preferred = freq_source;
end
