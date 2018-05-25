function h = vertical_in_ylim(t, varargin)
%HEM.PLOT.VERTICAL_IN_YLIM Plot a vertical line within axis dimensions.
%   This function will plot a vertical line within the current axis's (via
%   GCA) ylim. This is handy to mark the time of an event, such as LPS
%   injection or hormone infusion.
%
%   HEM.PLOT.VERTICAL_IN_YLIM(T) plots a vertical line at T. This line can
%   be styled similar to how a line made using PLOT can be styled by
%   passing any further arguments after T.
%
%   See also PLOT, GCA.

yl = get(gca, 'ylim');
hold on;
h = plot([t, t], yl, varargin{:});
ylim(yl);