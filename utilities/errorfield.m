function [varargout] = errorfield(x,y,err,varargin)
%plot error field: 
%h = errorfield(x,y,err,linespec,property name, property value...)
%plots one vector

ttl = [];
xlab = [];
ylab = [];

if nargin > 3
    linespec = varargin{1};
else
    linespec = 'b-';
end

if nargin > 4
    for a = 2:2:length(varargin)
        switch(varargin{a})
            case 'title'
                ttl = varargin{a+1};
            case 'xlabel'
                xlab = varargin{a+1};
            case 'ylabel'
                ylab = varargin{a+1};
            otherwise
                prop = varargin{a};
                propval = varargin{a+1};
        end
    end
end

if size(x,1) > size(x,2)
    x = x';
end
if size(y,1) > size(y,2)
    y = y';
end
if size(err,1) > size(err,2)
    err = err';
end

if (size(x,1) ~= 1) || (size(y,1) ~= 1)
    disp('Error: must be vector input');
    return;
end


h.h2 = fill([x fliplr(x)],[y+err fliplr(y-err)],linespec(1),'FaceAlpha',0.3,'EdgeColor',linespec(1),'EdgeAlpha',0.5);
holdstate =  ishold;
hold on;
h.h1 = plot(x,y,linespec,'LineWidth',2);
%hold off;
if ishold == 0
    hold off;
end
title(ttl);
xlabel(xlab);
ylabel(ylab);

%nargout

if nargout > 0
    varargout{1} = h;
end