function save_mat_variable_atomic(mat_file, variable_name, value)
%SAVE_MAT_VARIABLE_ATOMIC Save one MAT variable through a temporary file.

if nargin < 1 || isempty(mat_file)
    error('mat_file is required.');
end
if nargin < 2 || isempty(variable_name)
    error('variable_name is required.');
end

mat_file = char(string(mat_file));
variable_name = char(string(variable_name));
[mat_dir, ~, ~] = fileparts(mat_file);
if exist(mat_dir, 'dir') ~= 7
    mkdir(mat_dir);
end

use_short_tmp = strlength(mat_file) > 240 || strlength(mat_dir) > 220;
if use_short_tmp
    tmp_dir = fullfile(tempdir, 'koopman_mat_shortpath');
    if exist(tmp_dir, 'dir') ~= 7
        mkdir(tmp_dir);
    end
else
    tmp_dir = mat_dir;
end

tmp_file = [tempname(tmp_dir), '.mat'];
cleanup_tmp = onCleanup(@() local_delete_if_exists(tmp_file));

payload = struct();
payload.(variable_name) = value;
save(tmp_file, '-struct', 'payload', '-v7.3');
local_remove_existing(mat_file);

if use_short_tmp
    [ok, msg] = copyfile(tmp_file, mat_file, 'f');
    err_id = 'save_mat_variable_atomic:CopyFailed';
else
    [ok, msg] = movefile(tmp_file, mat_file, 'f');
    err_id = 'save_mat_variable_atomic:MoveFailed';
end
if ~ok
    error(err_id, ...
        'Unable to replace MAT file after temporary save:\n%s\n%s', ...
        mat_file, msg);
end

clear cleanup_tmp
end


function local_remove_existing(mat_file)
if exist(mat_file, 'file') ~= 2
    return;
end

try
    delete(mat_file);
catch delete_err
    [mat_dir, ~, ext] = fileparts(mat_file);
    backup_file = [tempname(mat_dir), '_replaced', ext];
    [ok, msg] = movefile(mat_file, backup_file, 'f');
    if ~ok
        error('save_mat_variable_atomic:ExistingFileLocked', ...
            ['Unable to remove existing MAT file before replace:\n%s\n' ...
            'Delete error: %s\nMove error: %s'], ...
            mat_file, delete_err.message, msg);
    end
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
