function poincare(x, varargin)
%HEM.PLOT.POINCARE Plot a poincare map.
%   HEM.PLOT.POINCARE(X) plots vector X(2:end) versus vector X(1:end-1).
%
%   Any arguments after X (e.g setting the marker style) are passed
%   directly to PLOT. The default markers are black dots with a MarkerSize
%   of 12.
%
%   See also HEM.MODEL.IPFM, HEM.PLOT.POINCARE_ELLIPSE,
%   HEM.UTIL.SAVE_MOVIE_POINCARE, PLOT.

if nargin < 2
    varargin = {'k.', 'MarkerSize', 12};
end

x1 = x(1:end-1);
x2 = x(2:end);

plot(x1, x2, varargin{:});
xlabel('RR(i) (s)');
ylabel('RR(i+1) (s)');