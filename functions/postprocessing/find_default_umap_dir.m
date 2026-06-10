function umap_dir = find_default_umap_dir(repo_root)
%FIND_DEFAULT_UMAP_DIR Resolve the preferred local UMAP toolbox location.

umap_dir = '';
icpbr_root = fileparts(fileparts(fileparts(repo_root)));
candidate_dirs = { ...
    fullfile(icpbr_root, 'forNikos', 'umapFileExchange (4.4)', 'umap'), ...
    fullfile(icpbr_root, 'forNikos', 'net_fmri_tutorials', ...
    'umapFileExchange (4.4)', 'umap')};

for i = 1:numel(candidate_dirs)
    candidate = candidate_dirs{i};
    if exist(fullfile(candidate, 'run_umap.m'), 'file') == 2
        umap_dir = candidate;
        return;
    end
end
end
