function [R, source_event_file] = load_event_results(cfg, output_root, event_input)
% Load saved event-detection results from a struct or file.

if nargin < 2 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 3
    event_input = [];
end

source_event_file = '';

if isstruct(event_input)
    if isfield(event_input, 'DetectResults')
        R = event_input;
        source_event_file = '[struct input]';
        return;
    end
    error('event_input struct must contain DetectResults.');
end

if ischar(event_input) || isstring(event_input)
    source_event_file = char(event_input);
    S = load(source_event_file);
    if ~isfield(S, 'R')
        error('The event-result file does not contain variable R.');
    end
    R = S.R;
    return;
end

search_dir = fullfile(output_root, cfg.file_stem, 'event_detection');
pattern = fullfile(search_dir, [cfg.file_stem, '_bandpass_events_*.mat']);
L = dir(pattern);

if isempty(L)
    error('No event-result file matching %s was found.', pattern);
end

if numel(L) > 1
    is_short_name = false(numel(L), 1);
    for i = 1:numel(L)
        is_short_name(i) = ~isempty(regexp(L(i).name, ...
            ['^', regexptranslate('escape', cfg.file_stem), '_bandpass_events_[0-9]+bands\.mat$'], 'once'));
    end

    if sum(is_short_name) == 1
        L = L(is_short_name);
    else
        names = cell(numel(L), 1);
        for i = 1:numel(L)
            names{i} = fullfile(L(i).folder, L(i).name);
        end
        error('Multiple event-result files found. Please pass event_input explicitly:\n%s', strjoin(names, newline));
    end
end

source_event_file = fullfile(L(1).folder, L(1).name);
S = load(source_event_file);
if ~isfield(S, 'R')
    error('The event-result file does not contain variable R.');
end

R = S.R;
end
