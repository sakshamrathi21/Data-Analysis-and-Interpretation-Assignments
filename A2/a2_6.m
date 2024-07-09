% The following is the MATLAB code for Question 6 of Assignment 2 of CS 215
% course. The code has been commented extensively to facilitate the reader
% in getting deeper insights into what each part of the code does.

% The code begins from here:

I1 = double(imread('T1.jpg'));          % Reading image 1
I2 = double(imread('T2.jpg'));          % Reading image 2
I1_prime = 255 - double(imread('T1.jpg'));  % Negative of image 1

% Defining the range of shifts to be performed
t_x = 10;
t_x_range = -t_x:t_x;

% Calculating the maximum possible size of the array
max_possible_size = 2 * max(abs(t_x_range)) + 1;

% Initializing correlation_coefficients with zeros
correlation_coefficients = zeros(1, max_possible_size);
correlation_coefficients_prime = zeros(1, max_possible_size);

% Initializing QMI arrays
qmi_values = zeros(1, length(t_x_range));
qmi_values_prime = zeros(1, length(t_x_range));

% Looping over different shift values
for i = 1:length(t_x_range)
    tx = t_x_range(i);
    
    % Shifting the images
    shifted_image = shiftImage(I2, tx);
    shifted_image_prime = shiftImage(I1_prime, tx);
    
    % Calculate correlation coefficients using your custom function
    correlation_coefficients(i) = calculateCorrelationCoefficient(I1(:), shifted_image(:));
    correlation_coefficients_prime(i) = calculateCorrelationCoefficient(I1(:), shifted_image_prime(:));
    
    % Calculating QMI
    [rho, qmi] = compute_metrics(I1, shifted_image, tx);
    [rho_prime, qmi_prime] = compute_metrics(I1, shifted_image_prime, tx);
    
    qmi_values(i) = qmi;
    qmi_values_prime(i) = qmi_prime;
end

% Creating a figure with subplots
figure('Position', [200, 200, 1200, 900]); % The numbers are adjusted for better view.

% Plotting ρ vs. tx for I2
subplot(2, 2, 1);
plot(t_x_range, correlation_coefficients, 'o-');
xlabel('Shift Values (tx)');
ylabel('Correlation Coefficients (ρ)');
title('Correlation Coefficients vs. Shift Values (I2)');
grid on;

% Plotting QMI vs. tx for I2
subplot(2, 2, 2);
plot(t_x_range, qmi_values, 'o-');
xlabel('Shift Values (tx)');
ylabel('Quadratic Mutual Information (QMI)');
title('QMI vs. Shift Values (I2)');
grid on;

% Plotting ρ vs. tx for I2_prime
subplot(2, 2, 3);
plot(t_x_range, correlation_coefficients_prime, 'o-');
xlabel('Shift Values (tx)');
ylabel('Correlation Coefficients (ρ)');
title('Correlation Coefficients vs. Shift Values (I2 Prime)');
grid on;

% Plotting QMI vs. tx for I2_prime
subplot(2, 2, 4);
plot(t_x_range, qmi_values_prime, 'o-');
xlabel('Shift Values (tx)');
ylabel('Quadratic Mutual Information (QMI)');
title('QMI vs. Shift Values (I2 Prime)');
grid on;

sgtitle('Dependence Measures vs. Shift Values'); % Adding a main title to the figure window.

% Rest of the code remains the same

% Function to compute joint histogram, histograms, ρ, and QMI
function [rho, qmi] = compute_metrics(I1, I2, ~)
    bin_width = 10;
    bins = 0:bin_width:255;
    num_bins = numel(bins) - 1;
    
    % Initializing histograms
    hist_I1 = zeros(1, num_bins);
    hist_I2 = zeros(1, num_bins);
    joint_hist = zeros(num_bins, num_bins);
    
    % Computing histograms
    for x = 1:num_bins
        hist_I1(x) = sum(I1(:) >= bins(x) & I1(:) < bins(x+1));
        hist_I2(x) = sum(I2(:) >= bins(x) & I2(:) < bins(x+1));
        for y = 1:num_bins
            joint_hist(x, y) = sum((I1(:) >= bins(x) & I1(:) < bins(x+1)) & ...
                                   (I2(:) >= bins(y) & I2(:) < bins(y+1)));
        end
    end
    
    % Normalizing histograms
    hist_I1 = hist_I1 / sum(hist_I1);
    hist_I2 = hist_I2 / sum(hist_I2);
    joint_hist = joint_hist / sum(joint_hist(:));
    
    % Computing ρ (correlation coefficient)
    bin_centers = (bins(1:end-1) + bins(2:end)) / 2; 
    mean_I1 = sum(bin_centers .* hist_I1);
    mean_I2 = sum(bin_centers .* hist_I2);
    cov_I1I2 = sum(sum((bin_centers' - mean_I1) .* (bin_centers - mean_I2) .* joint_hist));
    var_I1 = sum((bin_centers - mean_I1).^2 .* hist_I1);
    var_I2 = sum((bin_centers - mean_I2).^2 .* hist_I2);
    rho = cov_I1I2 / sqrt(var_I1 * var_I2);
    
    % Computing QMI (quadratic mutual information)
    pI1I2 = joint_hist;
    pI1 = hist_I1;
    pI2 = hist_I2;
    qmi = sum(sum((pI1I2 - (pI1' * pI2)).^2));
end

function shiftedImage = shiftImage(image, shiftAmount)
    [rows, cols, channels] = size(image);
    shiftedImage = zeros(rows, cols, channels, class(image));

    for channel = 1:channels
        if shiftAmount > 0
            % Positive shift to the right
            shiftedImage(:, 1:shiftAmount, channel) = 0;
            shiftedImage(:, shiftAmount+1:cols, channel) = image(:, 1:cols-shiftAmount, channel);
            % shiftedImage(:, 1:cols-shiftAmount, channel) = image(:, shiftAmount+1:cols, channel);
            % shiftedImage(:, cols-shiftAmount+1:cols, channel) = 0; % Fill the rightmost columns with zeros
        elseif shiftAmount < 0
            % Negative shift to the left
            shiftedImage(:, cols+shiftAmount+1:cols, channel) = 0;
            shiftedImage(:, 1:cols+shiftAmount, channel) = image(:, 1-shiftAmount:cols, channel);
            % shiftedImage(:, -shiftAmount+1:cols, channel) = image(:, 1:cols+shiftAmount, channel);
            % shiftedImage(:, 1:-shiftAmount, channel) = 0; % Fill the leftmost columns with zeros
        else
            % No shift, copy the original image
            shiftedImage(:, :, channel) = image(:, :, channel);
        end
    end
end

function correlation_coefficient = calculateCorrelationCoefficient(image1, image2)
    % Calculate the means of the two images
    mean1 = mean(image1);
    mean2 = mean(image2);

    % Calculate the Pearson correlation coefficient
    numerator = sum((image1 - mean1) .* (image2 - mean2));
    denominator = sqrt(sum((image1 - mean1).^2) * sum((image2 - mean2).^2));
    correlation_coefficient = numerator / denominator;
end

