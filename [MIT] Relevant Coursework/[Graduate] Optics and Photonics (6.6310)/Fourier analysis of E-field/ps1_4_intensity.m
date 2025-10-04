% Constants
c = 3e8;  
lambda_0 = 850e-9;
T_FWHM = 10e-15;  
omega_0 = 2 * pi * c / lambda_0;  
tau = T_FWHM / (2 * sqrt(2 * log(2))); 
samples = 1e3;
L = 10e-3;
um_to_m = 1e-6;
m_to_um = 1e6;
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
C1 = (0.0684043)^2;
C2 = (0.1162414)^2;
C3 = (9.896161)^2;

% Time
% "appropriate scale" based on 1.3
t = linspace(-25e-15, 25e-15, samples);

% time step; thought it would be helpful for the ifft circle shifting.
dt = t(2) - t(1);

% Gaussian pulse in time domain
% Given formula
E_0_t = exp(-(t.^2) / (2 * tau^2)) .* exp(1j * omega_0 * t);

% Fourier transform to get the spectrum
%FFT puts these in f (not omega)
E_0_omega = fftshift(fft(E_0_t));

frequencies = linspace(-1e14, 1e14, samples) + c/lambda_0;

% Lambda
lambda_m = linspace(0.6, 1.6, samples) * um_to_m;
lambda_um = lambda_m * m_to_um;

% n
n = sqrt(1 + (B1 .* lambda_um.^2 ./ (lambda_um.^2 - C1)) + (B2 .* lambda_um.^2 ./ (lambda_um.^2 - C2)) + (B3 .* lambda_um.^2 ./ (lambda_um.^2 - C3)));

% Beta
beta = (2 * pi * n) ./ lambda_m;  % Beta(ω) = (ω/c) * n(ω)

% Apply beta to FFT of initial E-field
E_1_omega = E_0_omega .* exp(-1j * beta * L);

% Inverse Fourier transform to get the time-domain pulse after propagation
E_1_t_L = ifft(ifftshift(E_1_omega));

% -1.54e-14 arbitrarily determined from trying to get the pulse to be
% symmetric on the vertical axis
N_shift = round(-1.47e-14 / dt);
E_1_t_L_shifted = circshift(E_1_t_L, N_shift);

% E-fields for comparison based on the dispersion
figure;
subplot(3, 1, 1);
plot(t * 1e15, real(E_1_t_L_shifted), 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('E_1(t, L) Shifted');
title('Shifted Dispersed Electric Field in Time Domain');
grid on;

% Original Output E-field
subplot(3, 1, 2);
plot(t * 1e15, real(E_1_t_L), 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('E_1(t, L)');
title('Original Dispersed Electric Field in Time Domain');
grid on;

%Original Input E-field
subplot(3, 1, 3);
plot(t * 1e15, real(E_0_t), 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('E_0(t)');
title('Original Input Electric Field in Time Domain');
grid on;

% Intensities
figure;
subplot(2, 1, 1);
plot(t * 1e15, abs(E_1_t_L_shifted).^2, 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('Intensity (I/I_0)');
title('Intensity of Dispersed Pulse (Linear Scale)');
grid on;

subplot(2, 1, 2);
plot(t * 1e15, 10*log10(abs(E_1_t_L_shifted).^2), 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('Intensity (dB)');
title('Intensity of Dispersed Pulse (Logarithmic Scale)');
grid on