function save_bold_post_mat(post_file, BOLD_POST)
%SAVE_BOLD_POST_MAT Safely save a pipeline 7 BOLD_POST MAT file.

if nargin < 1 || isempty(post_file)
    error('post_file is required.');
end
post_file = char(string(post_file));
[post_dir, ~, ~] = fileparts(post_file);
if exist(post_dir, 'dir') ~= 7
    mkdir(post_dir);
end

tmp_file = [tempname(post_dir), '.mat'];
cleanup_tmp = onCleanup(@() local_delete_if_exists(tmp_file));

save(tmp_file, 'BOLD_POST', '-v7.3');
local_remove_existing(post_file);

[ok, msg] = movefile(tmp_file, post_file, 'f');
if ~ok
    error('save_bold_post_mat:MoveFailed', ...
        'Unable to replace BOLD_POST file after temporary save:\n%s\n%s', ...
        post_file, msg);
end

clear cleanup_tmp
end


function local_remove_existing(post_file)
if exist(post_file, 'file') ~= 2
    return;
end

try
    delete(post_file);
catch delete_err
    [post_dir, ~, ext] = fileparts(post_file);
    backup_file = [tempname(post_dir), '_replaced', ext];
    [ok, msg] = movefile(post_file, backup_file, 'f');
    if ~ok
        error('save_bold_post_mat:ExistingFileLocked', ...
            ['Unable to remove existing BOLD_POST file before replace:\n%s\n' ...
            'Delete error: %s\nMove error: %s'], ...
            post_file, delete_err.message, msg);
    end
end

if exist(post_file, 'file') == 2
    error('save_bold_post_mat:ExistingFileStillPresent', ...
        'Existing BOLD_POST file could not be removed before replace:\n%s', ...
        post_file);
end
end


function local_delete_if_exists(path_in)
if exist(path_in, 'file') == 2
    try
        delete(path_in);
    catch
    end
end
end
