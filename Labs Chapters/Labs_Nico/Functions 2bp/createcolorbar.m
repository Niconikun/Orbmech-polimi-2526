function createcolorbar(ax)
%CREATECOLORBAR Create a colorbar that always matches data limits
%   AX: handle to the axes whose colorbar will be adjusted.

    % Create colorbar for the given axes
    cb = colorbar(ax);

    % Get current color limits of the axes

    % Set the colorbar limits explicitly to match the data range
    cb.Limits = clim(ax);

    % Make sure colorbar updates automatically when data changes
    addlistener(ax, 'CLim', 'PostSet', @(src, evt) set(cb, 'Limits',ax.CLim));

end