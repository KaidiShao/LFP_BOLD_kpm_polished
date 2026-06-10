function path_out = save_bold_reskoopnet_postprocessing_figure(fig, fig_dir, name, params)
%SAVE_BOLD_RESKOOPNET_POSTPROCESSING_FIGURE Save one pipeline 7 figure artifact.

path_out = '';
if isempty(fig) || ~isvalid(fig)
    return;
end

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    try
        disableDefaultInteractivity(axs(i));
    catch
    end
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end

if params.save_png
    path_out = fullfile(fig_dir, name);
    local_save_png(fig, path_out, params);
end
if params.save_fig
    [~, stem] = fileparts(name);
    savefig(fig, fullfile(fig_dir, [stem, '.fig']));
end
end

function local_save_png(fig, path_out, params)
method = local_get_param(params, 'figure_export_method', 'print');
method = lower(string(method));
if isempty(method) || strlength(method) == 0
    method = "print";
end
if method == "exportgraphics"
    exportgraphics(fig, path_out, 'Resolution', params.resolution);
    return;
end

try
    renderer = local_get_param(params, 'figure_renderer', '');
    if ~isempty(renderer)
        set(fig, 'Renderer', char(renderer));
    end
    drawnow limitrate;
    print(fig, path_out, '-dpng', sprintf('-r%d', params.resolution));
catch
    exportgraphics(fig, path_out, 'Resolution', params.resolution);
end
end

function value = local_get_param(params, field_name, default_value)
value = default_value;
if isstruct(params) && isfield(params, field_name)
    value = params.(field_name);
end
end
