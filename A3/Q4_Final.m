% The following is the MATLAB code for Question 4 of Assignment 3 of CS 215
% course. The code has been commented extensively to facilitate the reader
% in getting deeper insights into what each part of the code does.

clc; clf; clear;
%rng('shuffle'); % to generate a new set of random values each time

% Generating independent samples from N(0, 16)
samples = 4*randn(1,1000);

% Creating disjoint subsets T and V of A
index_perm = randperm(1000);
T_size = 750; 
T = samples(index_perm(1:T_size));
V = samples(index_perm((T_size+1):1000));

% set of possible sigma values
sigma_values = [0.001, 0.1, 0.2, 0.9, 1, 2, 3, 5, 10, 20, 100];

% ll_val is the set Collection of Log of LL Values
ll_val = [];
for sigma = sigma_values
    temp=[];
    for v = V
        x = ((v-T).^2);
        al = log(sum(exp(-x/(2*sigma^2))));
        temp = [temp al];
    end
   
    ll_val = [ll_val sum(temp)-(1000-T_size)*log(T_size*sigma*sqrt(2*pi))];
end

% Plotting graph of LL vs log(sigma) 
plot(log(sigma_values), ll_val);
title('LL vs log(sigma)');
[val, pos] = max(ll_val);
fprintf('Best LL is given for sigma = %0.3f\n', sigma_values(pos));
grid on;

best_sigma = sigma_values(pos);

x_val = -8:0.1:8;

% pld is the set of values of the pdf for different values of x in x_val
pld = [];
for x = x_val
    a = ((x - T).^2);
    al = (sum(exp(-a/(2*best_sigma^2))))/(T_size*best_sigma*sqrt(2*pi));
    pld = [pld al];
end

figure;
plot(x_val, pld, 'r-', 'LineWidth', 2);
hold on;
true_density = normpdf(x_val, 0, 4); % because values of T were drawn from N(0, 4)
plot(x_val ,true_density, 'b--', 'LineWidth', 2);

xlabel('Values');
ylabel('Probability Density');
title('Estimated Density vs. True Density');
legend('Estimated Density', 'True Density');
grid on;

ll_val=[];
for sigma=sigma_values
    temp = [];
    for v=V
        ph = sum(exp(-((v-T).^2)/(2*sigma^2)))/(T_size*sigma*sqrt(2*pi));
        x = (exp(-v^2/32)/(4*sqrt(2*pi)) - ph)^2;
        temp = [temp x];
    end
    ll_val = [ll_val sum(temp)];
end

figure;
plot(log(sigma_values), ll_val);
title('D vs log(sigma)');
[v,po] = min(ll_val);
fprintf('Best D is given for sigma = %.3f\n', sigma_values(po));
fprintf('D value for sigma giving best LL = %.5f\n', ll_val(pos));
grid on;

D_best_sigma = sigma_values(po);

pldd = [];
for x = x_val
    a = ((x - T).^2);
    al = (sum(exp(-a/(2*D_best_sigma^2))))/(T_size*D_best_sigma*sqrt(2*pi));
    pldd = [pldd al];
end

figure;
plot(x_val, pldd, 'r-', 'LineWidth', 2);
hold on;
true_density = normpdf(x_val, 0, 4); % because values of T were drawn from N(0, 4)
plot(x_val ,true_density, 'b--', 'LineWidth', 2);

xlabel('Values');
ylabel('Probability Density');
title('Estimated Density vs. True Density');
legend('Estimated Density', 'True Density');
grid on;