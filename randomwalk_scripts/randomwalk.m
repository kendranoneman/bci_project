clc; clearvars; 
fontSize = 20;
format compact;

numberOfSteps = 75;
mu = 3;
sigma = 1; 

% Get all the random step directions in advance (for speed) rather than in the loop.
deltax = ones(1, numberOfSteps);
%deltay = rand(1, numberOfSteps) - 0.5;
deltay = normrnd(mu, sigma,[1, numberOfSteps]);
% Define starting location.
x = zeros(1, numberOfSteps);
y = zeros(1, numberOfSteps);
text(x, y, '1', 'fontSize', fontSize);
hold on;
% Now do the steps in a random direction.
for step = 2 : numberOfSteps
	% Walk in the x direction.
	x(step) = x(step-1) + deltax(step);
	% Walk in the y direction.
	y(step) = y(step-1) + deltay(step);
	% Now plot the walk so far.
	plot(x(1:step), y(1:step), 'bo-', 'LineWidth', 2);
	textLabel = sprintf('%d', step);
	text(x(step), y(step), textLabel, 'fontSize', fontSize);
end

% Plot x and y axes
line(xlim, [0,0], 'Color', 'k', 'LineWidth', 2);
line([0,0], ylim, 'Color', 'k', 'LineWidth', 2);

% Mark the first point in red.
hold on;
plot(x(1), y(1), 'rs', 'LineWidth', 2, 'MarkerSize', 25);
textLabel = '1';
text(x(1), y(1), textLabel, 'fontSize', fontSize);
grid on;

% Mark the last point in red.
plot(x(end), y(end), 'rs', 'LineWidth', 2, 'MarkerSize', 25);
caption = sprintf('Random Walks with %d steps', numberOfSteps);
title(caption, 'FontSize', fontSize);
xlabel('Time from BCI start (1 bin = 50 ms)', 'FontSize', fontSize);
ylabel('Avg. Distance', 'FontSize', fontSize);

% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0, 0.04, 1, 0.96]);

% Calculate the distance from the origin.
distanceFromOrigin = hypot(x(end), y(end));
message = sprintf('Done with demo!\nDistance of endpoint from origin = %.3f', distanceFromOrigin);
msgbox(message);
