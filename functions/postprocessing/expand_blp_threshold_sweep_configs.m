function [cfgs, ratio_tags] = expand_blp_threshold_sweep_configs(cfg)
%EXPAND_BLP_THRESHOLD_SWEEP_CONFIGS Clone a threshold config over ratios.

if isfield(cfg, 'threshold_ratio_sweep') && ~isempty(cfg.threshold_ratio_sweep)
    ratios = double(cfg.threshold_ratio_sweep(:)).';
else
    ratios = double(cfg.threshold_ratio);
end

if isempty(ratios) || any(~isfinite(ratios))
    error('threshold_ratio_sweep must contain finite numeric values.');
end

if ~isfield(cfg, 'threshold_ratio_sweep')
    cfg.threshold_ratio_sweep = [];
end
if ~isfield(cfg, 'threshold_sweep_tag')
    cfg.threshold_sweep_tag = '';
end

cfgs = repmat(cfg, 1, numel(ratios));
ratio_tags = strings(1, numel(ratios));
use_subdirs = numel(ratios) > 1;

for i = 1:numel(ratios)
    ratio = ratios(i);
    ratio_tag = local_threshold_ratio_tag(ratio, cfg);
    ratio_tags(i) = string(ratio_tag);

    cfg_i = cfg;
    cfg_i.threshold_ratio = ratio;
    cfg_i.threshold_ratio_sweep = [];
    cfg_i.threshold_sweep_tag = ratio_tag;

    if use_subdirs
        if isfield(cfg_i, 'save_dir') && ~isempty(cfg_i.save_dir)
            cfg_i.save_dir = fullfile(cfg_i.save_dir, ratio_tag);
        end
        if isfield(cfg_i, 'figure_dir') && ~isempty(cfg_i.figure_dir)
            cfg_i.figure_dir = fullfile(cfg_i.figure_dir, ratio_tag);
        end
    end

    cfgs(i) = cfg_i;
end
end


function tag = local_threshold_ratio_tag(ratio, cfg)
mode = 'ratio';
if isfield(cfg, 'threshold_mode') && ~isempty(cfg.threshold_mode)
    mode = lower(char(cfg.threshold_mode));
end

switch mode
    case 'quantile'
        prefix = 'q';
    case 'maxfrac'
        prefix = 'maxfrac';
    case 'meanplusstd'
        prefix = 'mps';
    otherwise
        prefix = 'ratio';
end

if abs(ratio) < 100
    value = sprintf('%03d', round(ratio * 100));
else
    value = regexprep(sprintf('%.6g', ratio), '[^\w]+', 'p');
end

tag = sprintf('%s%s', prefix, value);
end
