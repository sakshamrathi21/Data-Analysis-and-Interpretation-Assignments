clear; clc; clf;
%rng(42);

x = -3.0:0.02:3.0;
y = 6.5 * sin(2.1 * x + pi/3);

xlabel('x');
ylabel('y');
plot(x, y, 'DisplayName', 'Orignal y', color='black');

% Randomly choose a fraction f from y
f = 0.3;  
len_y = length(y);
rand_indices = randperm(len_y, floor(len_y*f));

z = y;

for i=1:length(rand_indices)
    corruption = 100 + 20*rand;
    z(rand_indices(i)) = y(rand_indices(i)) + corruption;
end

%%% Filtering corrupted subset of y %%%

% Moving median filtering
y_median = filter(@median, z);

hold on;
plot(x, y_median, 'DisplayName', 'y median', color='red');

% Moving average filtering
y_mean = filter(@mean, z);

hold on;
plot(x, y_mean, 'DisplayName', 'y mean', color='green');

% Moving quartile filtering (1st quartile = 25th percentile)
y_quantile = filter(@quartile, z);

hold on;
plot(x, y_quantile, 'DisplayName', 'y quartile', color='cyan');

legend

% Computing relative mean squared error between y_filtered and y
disp("rms median   = " + rms(y_median, y));
disp("rms mean     = " + rms(y_mean, y));
disp("rms quantile = " + rms(y_quantile, y));

% Helper functions
function q = quartile(A)
    q = quantile(A, 0.25);
end

function result = filter(central_tendency, z)
    len_z = length(z);
    y_central_tendency = zeros(1, len_z);

    for i=1:len_z
        len_N = min(i-1, 8) + min(8, len_z-i+1) + 1;
        N = zeros(1, len_N);
        
        for j=1:len_N
            N(j) = z(j + i-min(i, 10)); 
        end
    
        y_central_tendency(i) = central_tendency(N);
    end

    result = y_central_tendency;
end

function result = rms(y_hat, y)
    diff = y_hat - y;
    result = dot(diff, diff) / dot(y, y);
end

