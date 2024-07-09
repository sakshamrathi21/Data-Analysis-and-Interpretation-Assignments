N_MAX = 100000;
n_values = [10, 100, 1000, 10000, N_MAX]; % Variety of n values

figure('Position', [100, 100, 1200, 800]);

for i = 1:length(n_values)
    n = n_values(i);
    array_expectations = zeros(1, n);

    for j = 1:n
        array_expectations(j) = E_X_N(j);
    end

    % Create a subplot grid
    subplot(2, 3, i); % 2x3 grid of subplots

    % Plot the values vs. n
    plot(1:n, array_expectations);
    xlabel('n');
    ylabel('E(X^{(n)})');
    title(['Expectation of X^{(n)} vs. n (n = ' num2str(n) ')']);
    grid on;
end

% Generate some sample data for n and E(XN)
n_values = 1:N_MAX; % Adjust the range as needed
EXN = zeros(1, N_MAX); % Initialize the array
for n = n_values
    EXN(n) = E_X_N(n);
end

% % Calculate the other two curves
% nlogn_over_4 = (n_values .* log(n_values)) / 4;
% twonlogn = 2 * (n_values .* log(n_values));
% 
% % Create a subplot
% subplot(2, 3, 6);
% 
% % Plot the three curves
% plot(n_values, EXN, 'b', 'LineWidth', 2);
% hold on;
% plot(n_values, nlogn_over_4, 'g', 'LineWidth', 2);
% plot(n_values, twonlogn, 'r', 'LineWidth', 2);
% 
% % Add labels and legend
% xlabel('n');
% ylabel('Value');
% title('E(XN), nlogn/4, and 2nlogn Curves');
% legend('E(XN)', 'nlogn/4', '2nlogn');
% 
% % Hold off to stop overlaying plots
% hold off;

% Adjust the layout
sgtitle('Variety of Expectations vs. n');

function [expectation] = E_X_N(n)
    expectation = 0;
    for i = 1:n
        expectation = expectation + 1/i;
    end
    expectation = expectation * n;
end
