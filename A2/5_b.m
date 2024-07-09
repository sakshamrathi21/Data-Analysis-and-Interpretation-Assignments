clc; clf; clear;

values = 1:5;
probabs = [0.05, 0.4, 0.15, 0.3, 0.1];

Ns = [5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000];

% Set the desired width and height in pixels
width = 600;
height = 300;

% Create a new figure with the specified width and height
fig = figure('Position', [100, 100, width, height]);

for N = Ns 
    % Generate the Xs data for the current N
    Xs = zeros(N, 1);
    for i = 1:N
        X = randsample(values, N, true, probabs);
        X_avg = mean(X);
        Xs(i) = X_avg;
    end
    
    % Create a subplot for each N
    subplot(2, 5, find(Ns == N)); % 2 rows, 5 columns

    % Plot the ECDF for X_avg's
    ecdf(Xs);

    mu = mean(Xs);
    sigma = std(Xs);
    
    x = linspace(mu - 4*sigma, mu + 4*sigma, 1000);
    p = normcdf(x, mu, sigma);
    
    % Plot the Normal CDF
    hold on;
    plot(x, p);
    
    % Set the x-axis limits
    xlim([mu - 4*sigma, mu + 4*sigma]);

    % labels
    title(['N = ', num2str(N)]);
    
    if find(Ns == N) > 5
        xlabel('X');
    end
    if mod(find(Ns == N), 5) == 1
        ylabel('CDF');
    end
end

% Add the custom legend entry for Normal CDF to the last subplot
subplot(2, 5, 10);
legend('Normal CDF', 'Location', 'best');

sgtitle('Empirical vs. Normal CDF');
