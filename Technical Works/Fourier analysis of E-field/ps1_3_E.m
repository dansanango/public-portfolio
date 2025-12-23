% Constants
c = 3e8;
lambda_0 = 850e-9;
T_FWHM = 10e-15;
omega_0 = 2 * pi * c / lambda_0;
tau = T_FWHM / (2 * sqrt(2 * log(2)));
samples = 1e7;

% Time
t = linspace(-25e-15, 25e-15, samples);

% Gaussian pulse in time domain
E_0_t = exp(-(t.^2) / (2 * tau^2)) .* exp(1j * omega_0 * t);

% Intensity of the pulse
% Technically this should be multiplied by 1/(2*eta), but the Pset later
% says I~E^2 so I'm going to assume this variable is essentially negligible
% for the purpose of this question. This will go for all later intensity
% calculations.
I_t = abs(E_0_t).^2;

% Fourier transform to get the spectrum and bandwidth
% Note: This value is in frequency, not omega
E_0_omega = fftshift(fft(E_0_t));

%Covers frequency range {-1/2, 1/2}*1/time interval [t(2) - t(1)]
frequencies = linspace(-1e14, 1e14, samples) + c/lambda_0;

% Input E
figure;
subplot(3, 1, 1);
plot(t * 1e15, real(E_0_t), 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('E_0(t) (V/m)');
title('Input Electric Field');
grid on;

% Initial pulse intensity
subplot(3, 1, 2);
plot(t * 1e15, I_t, 'LineWidth', 1.5);
xlabel('Time (fs)');
ylabel('Intensity (I/I_0)');
title('Input Pulse Intensity');
grid on;

% Magnitude of the Fourier transform
subplot(3, 1, 3);
plot(frequencies * 1e-12, abs(E_0_omega), 'LineWidth', 1.5);
xlabel('Frequency (THz)');
ylabel('|E_0(omega)|');
title('Pulse in Frequency Domain');
xlim([(c/lambda_0) * 1e-12 - 1e-3, (c/lambda_0) * 1e-12 + 1e-3]);
grid on;
