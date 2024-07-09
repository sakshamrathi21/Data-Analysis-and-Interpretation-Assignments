clc; clf; clear;

values = 1:5;
probabs = [0.05, 0.4, 0.15, 0.3, 0.1];

Ns = [5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000];

figure('Name', 'Histograms for Different N');

for idx = 1:length(Ns)
    N = Ns(idx);
    
    % Generate the Xs data for the current N
    Xs = zeros(N, 1);
    for i = 1:N
        X = randsample(values, N, true, probabs);
        X_avg = mean(X);
        Xs(i) = X_avg;
    end

    % Create a subplot for each N
    subplot(2, 5, idx); % 2 rows, 5 columns, and select the current index
    histogram(Xs, 50);
    title(['N = ', num2str(N)]);
end

sgtitle('Histograms for Different N'); 