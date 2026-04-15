function results_root = get_project_results_root(repo_root)
%GET_PROJECT_RESULTS_ROOT Resolve the default writable results root.
%
%   results_root = GET_PROJECT_RESULTS_ROOT(repo_root)
%
% Priority:
%   1. Environment variable LFP_BOLD_KPM_RESULTS_ROOT
%   2. Windows default external drive path
%   3. repo_root/results fallback

if nargin < 1 || isempty(repo_root)
    this_file_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(this_file_dir));
end

env_root = strtrim(getenv('LFP_BOLD_KPM_RESULTS_ROOT'));
if ~isempty(env_root)
    results_root = env_root;
    return;
end

if ispc
    results_root = 'E:\LFP_BOLD_kpm_polished_results';
else
    results_root = fullfile(repo_root, 'results');
end
end
