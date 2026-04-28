function [C, source_consensus_file] = load_consensus_state_results(cfg, output_root, consensus_input)
% Load saved consensus-state results from a struct or file.

if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 3
    consensus_input = [];
end

if isstruct(consensus_input)
    if isfield(consensus_input, 'state_windows') && isfield(consensus_input, 'state_catalog')
        C = consensus_input;
        source_consensus_file = '[struct input]';
        return;
    end
    error('consensus_input struct must contain state_windows and state_catalog.');
end

if ischar(consensus_input) || isstring(consensus_input)
    source_consensus_file = char(consensus_input);
    S = load(source_consensus_file);
    if ~isfield(S, 'C')
        error('The consensus-state file does not contain variable C.');
    end
    C = S.C;
    return;
end

search_dir = fullfile(output_root, cfg.file_stem, 'consensus_states');
pattern = fullfile(search_dir, [cfg.file_stem, '_consensus_states_*.mat']);
L = dir(pattern);

if isempty(L)
    error('No consensus-state file matching %s was found.', pattern);
end

if numel(L) > 1
    names = cell(numel(L), 1);
    for i = 1:numel(L)
        names{i} = fullfile(L(i).folder, L(i).name);
    end
    error('Multiple consensus-state files found. Please pass consensus_input explicitly:\n%s', strjoin(names, newline));
end

source_consensus_file = fullfile(L(1).folder, L(1).name);
S = load(source_consensus_file);
if ~isfield(S, 'C')
    error('The consensus-state file does not contain variable C.');
end

C = S.C;
end
