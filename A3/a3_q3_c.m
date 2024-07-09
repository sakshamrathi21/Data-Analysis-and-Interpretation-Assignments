filename = 'XYZ.txt';
delimiter = ',';
data = dlmread(filename, delimiter);

X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

P11 = sum(X.^2);
P12  = sum(X.*Y);
P13 = sum(X);
P21 = P12;
P22  = sum(Y.^2);
P23 = sum(Y);
P31 = P13;
P32 = P23;
P33  = 2000;
P = [P11, P12, P13; P21, P22, P23; P31, P32, P33];

q1 = sum(X.*Z);
q2 = sum(Z.*Y);
q3  =sum(Z);
q = [q1; q2; q3];

result = P\q;

noise = Z - (result(1) * X + result(2) * Y + result(3));
noise_mean = sum(noise) / 2000; 
noise_variance = sum((noise - noise_mean).^2) / (1999); 

text = "The predicted equation of the plane is " + result(1) + "x + " + result(2) + "y + " + result(3);
disp(text);

text1 = "The predicted noise variance is " + noise_variance;
disp(text1);

% Define the data points or equation of the plane
[X, Y] = meshgrid(-10:0.1:10, -10:0.1:10);
Z = result(1)*X + result(2)*Y + result(3); % Example equation of the plane: Z = 2*X + 3*Y + 5

% Create the 3D surface plot
figure;
surf(X, Y, Z);

% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Plane Plot');

% Adjust the view for better visualization
view(30, 30);
