function processed_root = get_project_processed_root()
%GET_PROJECT_PROCESSED_ROOT Resolve the default processed-data root.
%
% Priority:
%   1. Environment variable LFP_BOLD_KPM_PROCESSED_ROOT
%   2. Windows default external drive path
%   3. POSIX default mounted external drive path

env_root = strtrim(getenv('LFP_BOLD_KPM_PROCESSED_ROOT'));
if ~isempty(env_root)
    processed_root = env_root;
    return;
end

if ispc
    processed_root = 'E:\DataPons_processed\';
else
    processed_root = '/mnt/e/DataPons_processed/';
end
end
