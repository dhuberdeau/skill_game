function varargout = plotPath_fine(paths)
% function plotPath(paths)
%
% Plot the path on the current axes (or new if none are active)
tParm = 0:.001:1;
wholePath = nan(6, length(paths.Paths)*length(tParm));
lenPath = nan(3,length(paths.Paths));
iwp = 1;
for i = 1:length(paths.Paths)
    l = renderCurve([paths.Paths(i).c_Left(:,1), paths.Paths(i).c_Left(:,3)]', tParm);
    r = renderCurve([paths.Paths(i).c_Right(:,1), paths.Paths(i).c_Right(:,3)]', tParm);
    c = renderCurve([paths.Paths(i).c_Center(:,1), paths.Paths(i).c_Center(:,3)]', tParm);
    wholePath(1:2,iwp:(iwp+length(tParm)-1)) = l;
    wholePath(3:4,iwp:(iwp+length(tParm)-1)) = r;
    wholePath(5:6,iwp:(iwp+length(tParm)-1)) = c;
    lenPath(1,i) = sum(sqrt(diff(l(1,:)).^2 + diff(l(2,:)).^2));
    lenPath(2,i) = sum(sqrt(diff(r(1,:)).^2 + diff(r(2,:)).^2));
    lenPath(3,i) = sum(sqrt(diff(c(1,:)).^2 + diff(c(2,:)).^2));
    iwp = iwp+length(tParm);
    plot(l(1,:), l(2,:), 'r'); hold on;
    plot(c(1,:), c(2,:), 'Color', [.5 .5 .5]);
    plot(r(1,:), r(2,:), 'g');
    plot(r(1,end), r(2,end), 'go');
    text(l(1,1), l(2,1), num2str(i));
end
varargout{1} = wholePath;
varargout{2} = lenPath;