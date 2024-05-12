function [h] = plot_pose(T, color)


h = cell(4,1);

t = T(1:3,4);
R = T(1:3, 1:3);

% plot 2d x y z
hold on
h{1} = plot3( t(1), t(2), t(3), '.','Color', color, 'MarkerSize', 1); hold on;

%draw axis
axis = 5;

xdir = t + axis*R(:, 1);      %This is just the point of the local x-axis times 5
xdir = [t, xdir];
ydir = t + axis*R(:, 2);      % plot is just plot the point of local position, and the point of local x-axis times 5
ydir = [t, ydir];
zdir = t + axis*R(:, 3);
zdir = [t, zdir];

h{2} = plot3(xdir(1,:), xdir(2,:), xdir(3,:), 'Color', 'red', 'LineWidth',2);  hold on;
h{3} = plot3(ydir(1,:), ydir(2,:), ydir(3,:), 'Color', 'green','LineWidth',2); hold on;
h{4} = plot3(zdir(1,:), zdir(2,:), zdir(3,:), 'Color', 'blue', 'LineWidth',2); hold on;


end