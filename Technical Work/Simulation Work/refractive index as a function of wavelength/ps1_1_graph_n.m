% Constants
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
C1 = (0.0684043)^2;
C2 = (0.1162414)^2;
C3 = (9.896161)^2;
samples = 1e3;

% Wavelength from 0.6[um] to 1.6[um]
% the equation is given to us with lambda in terms of microns, so no
% scaling necessary
lambda = linspace(0.6, 1.6, samples);

% calculating n
n_squared = 1 + (B1*lambda.^2 ./ (lambda.^2 - C1)) + (B2*lambda.^2 ./ (lambda.^2 - C2)) + (B3*lambda.^2 ./ (lambda.^2 - C3));
n = sqrt(n_squared);

% Graphing
figure;
plot(lambda, n, 'LineWidth', 2);
xlabel('Wavelength (microns)');
ylabel('Refractive Index n(\lambda)');
title('Refractive Index vs Wavelength for Fused Silica');
grid on;
