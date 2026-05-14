function S = load_mat_file_with_short_path(file, varargin)
%LOAD_MAT_FILE_WITH_SHORT_PATH Load MAT files, falling back through a short temp path.

file = char(string(file));
if strlength(file) > 240 && exist(file, 'file') == 2
    S = local_load_via_short_copy(file, varargin{:});
    return;
end

try
    S = load(file, varargin{:});
    return;
catch first_err
    if exist(file, 'file') ~= 2
        rethrow(first_err);
    end
end

S = local_load_via_short_copy(file, varargin{:});
end


function S = local_load_via_short_copy(file, varargin)
tmp_dir = fullfile(tempdir, 'koopman_mat_shortpath');
if exist(tmp_dir, 'dir') ~= 7
    mkdir(tmp_dir);
end
[~, ~, ext] = fileparts(file);
if isempty(ext)
    ext = '.mat';
end
tmp_file = [tempname(tmp_dir), ext];
cleanup_tmp = onCleanup(@() local_delete_if_exists(tmp_file));

[ok, msg] = copyfile(file, tmp_file, 'f');
if ~ok
    error('load_mat_file_with_short_path:CopyFailed', ...
        'Unable to copy MAT file to a short temporary path before load:\n%s\n%s', ...
        file, msg);
end
S = load(tmp_file, varargin{:});

clear cleanup_tmp
end


function local_delete_if_exists(path_in)
if exist(path_in, 'file') == 2
    try
        delete(path_in);
    catch
    end
end
end
