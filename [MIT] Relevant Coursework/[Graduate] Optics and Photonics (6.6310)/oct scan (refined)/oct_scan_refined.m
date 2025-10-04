% Objective of question: in the previous question, we used an IFT to get the scan results. However, we did not account
% for phase correction. Also, because wavenumber and wavelength are inversely proportional, we will have large clumps of
% data in some places and more sparse data in others. To account for this, interpolate wavelengths.

% Load data
data = load('BScanFringe 2.mat');

% Get table from loaded data
k_space_data = data.fringe_interpolated;

% Define parameter for second-order dispersion correction (adjust as needed)
% guessed a until i got a higher-quality image
a = -1.3e1;

% Get dimensions of k-space data
[num_k, num_x] = size(k_space_data);
k_axis = linspace(-pi, pi, num_k);  % Adjust range as needed for your data

% Compute the second-order phase correction term
phase_correction = exp(1j * a * k_axis.^2);
phase_correction = repmat(phase_correction.', 1, num_x);

% Apply the dispersion correction to k-space data
k_space_data_corrected = k_space_data .* phase_correction;

% Perform ifft and ifftshift to get z-space data
z_space_data_corrected = ifftshift(ifft(k_space_data_corrected, [], 1));
z_space_data_corrected = ifftshift(z_space_data_corrected, 2);

% Convert to magnitude, log scale, and thresholding as before
z_space_magnitude_corrected = abs(z_space_data_corrected);
log_z_space_corrected = log10(z_space_magnitude_corrected + 1e-3);

% Apply thresholds to log-scaled data
log_z_space_clipped_corrected = log_z_space_corrected;
log_z_space_clipped_corrected(log_z_space_corrected < lower_threshold) = lower_threshold;
log_z_space_clipped_corrected(log_z_space_corrected > upper_threshold) = upper_threshold;

figure;
imagesc(log_z_space_clipped_corrected);  % Corrected image
colormap(gray);
axis image;
title('Corrected B-Scan');
xlabel('x-pixel');
ylabel('y-pixel');
