% Define x-axis range (time from 0 to 20 seconds)
x = [0 20]; 

% Define y-axis value (constant flow rate)
y = [2108.813609 2108.813609];

% Plot the line
plot(x, y, 'b-', 'LineWidth', 1);

% Labels and title
xlabel('Time (s)');
ylabel('Flow Rate (cm^3/s)');
title('Constant Flow Rate Over Time');

% Set axis limits
xlim([0 20]);
ylim([2000 2200]); % Adjust if needed

% Grid for better visualization
grid on;
