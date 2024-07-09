clc; clf; clear;

values = 1:5;
probabs = [0.05, 0.4, 0.15, 0.3, 0.1];

Ns = [5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000];

MADs = zeros(1, length(Ns));
for i = 1:length(Ns)
    N = Ns(i);

    % Generate the Xs data for the current N
    Xs = zeros(N, 1);
    for j = 1:N
        X = randsample(values, N, true, probabs);
        X_avg = mean(X);
        Xs(j) = X_avg;
    end
    
    [f, x] = ecdf(Xs);

    mu = mean(Xs);
    sigma = std(Xs);
    
    p = normcdf(x, mu, sigma);
    
    MADs(i) = max(abs(p - f));
end

plot(log10(Ns), MADs);
title('Maximum Absolute Deviation (MAD) vs. Sample Size (N)');
xlabel('Sample Size (N)');
ylabel('MAD');
