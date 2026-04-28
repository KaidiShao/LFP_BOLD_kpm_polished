function [session_ids, session_channels] = collect_included_session_entries(cfg, default_channels)
%COLLECT_INCLUDED_SESSION_ENTRIES Collect included session IDs and optional channel selections.
%
% With one input:
%   session_ids = collect_included_session_entries(cfg)
%
% With two inputs:
%   [session_ids, session_channels] = collect_included_session_entries(cfg, default_channels)
%
% When channel output is requested, each returned session ID is paired with
% one channel-selection entry. Empty cfg.sessions(i).selected_channels falls
% back to default_channels.

session_ids = [];
session_channels = {};

use_channels = (nargout >= 2) || (nargin >= 2);
if use_channels && nargin < 2
    default_channels = [];
end

for i = 1:numel(cfg.sessions)
    if ~cfg.sessions(i).include
        continue;
    end

    this_session_ids = cfg.sessions(i).session_id(:);

    if use_channels
        if isfield(cfg.sessions(i), 'selected_channels') && ...
                ~isempty(cfg.sessions(i).selected_channels)
            this_channels = cfg.sessions(i).selected_channels;
        else
            this_channels = default_channels;
        end
    end

    for j = 1:numel(this_session_ids)
        session_ids(end+1, 1) = this_session_ids(j); %#ok<AGROW>
        if use_channels
            session_channels{end+1, 1} = this_channels; %#ok<AGROW>
        end
    end
end
