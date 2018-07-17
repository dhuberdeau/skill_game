function varargout = renderCurve(P, t, varargin)
% function f = renderCurve(P, t, varargin)
%
% Renders the 2-D bezier curve corresponding to the control points in P and
% the parameter t.
%
% P = [p1, p2, p3, ..., pk] (p_i = column vector)
%
% t = increasing values (e.g. 0:.01:1)
%
% variable inputs:
%   (1) H - a figure or axis handle
%   (2) plotStyle - a designation for the color and line style.
%    * must be in this order (i.e. handle then plot style
%
% David Huberdeau, 6/14/2012

%% open the figure
switch nargin
    case 0
        error('Not enough inputs');
    case 1
        error('Not enough inputs');
    case 2
        H = -1;
        plotColor = 'b-';
    case 3
        H = varargin{1};
        plotColor = 'b-';
    case 4
        H = varargin{1};
        plotColor = varargin{2};
    otherwise
        error('Too many inputs');
end

%% compute the curve's points from de Castejian formula:
f = nan(2, length(t));
n = size(P,2);
for c = 1:length(t)
    Pset = P;
    for i = 1:(n-1)
        Ptemp = Pset(:, 1:(end - 1));
        for j = 1:size(Ptemp,2)
            Ptemp(:,j) = (1-t(c))*Pset(:,j) + t(c)*Pset(:,j+1);
        end
        Pset = Ptemp;
    end
    f(:,c) = Pset(:,end);
end

%% render on graph
if H > 0
    figure(H); hold on;
    try
        plot(f(1,:), f(2,:), plotColor);
    catch errMsg
        warning(errMsg);
        plot(f(1,:), f(2,:));
    end
end

varargout{1} = f;