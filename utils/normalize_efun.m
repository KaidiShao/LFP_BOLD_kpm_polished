function efuns_norm = normalize_efun(efuns, mode)
% Normalize each eigenfunction column by its maximum absolute magnitude.
%
% mode:
%   'real' -> normalize real(efuns)
%   'abs'  -> normalize abs(efuns)

efuns_norm = zeros(size(efuns, 1), size(efuns, 2));

if strcmpi(mode, 'real')
    for n_efun = 1:size(efuns, 2)
        x = real(efuns(:, n_efun));
        denom = max(abs(efuns(:, n_efun)));
        if denom > 0
            efuns_norm(:, n_efun) = x ./ denom;
        else
            efuns_norm(:, n_efun) = x;
        end
    end

elseif strcmpi(mode, 'abs')
    for n_efun = 1:size(efuns, 2)
        x = abs(efuns(:, n_efun));
        denom = max(abs(efuns(:, n_efun)));
        if denom > 0
            efuns_norm(:, n_efun) = x ./ denom;
        else
            efuns_norm(:, n_efun) = x;
        end
    end
else
    error('Unknown mode %s. Use ''real'' or ''abs''.', mode);
end
