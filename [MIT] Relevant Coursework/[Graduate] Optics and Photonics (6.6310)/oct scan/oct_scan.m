% Objective of question: From raw eye scan data, generate an image of the
% eye. This code goes through the following logic:

% 1. an IFT of the output data
% 2. making the data easier to see (taking the log of the output, filtering
% values that are too large/small...

% load data
data = load('BScanFringe 2.mat');

% get table from loaded data
k_space_data = data.fringe_interpolated;

% ifft, then ifftshift so the y-axis looks proper
z_space_data = ifftshift(ifft(k_space_data, [], 1));

% ifftshift so the x-axis looks proper
z_space_data = ifftshift(z_space_data, 2);

% take magnitude
z_space_magnitude = abs(z_space_data);
% log, add 1e-3 to avoid log(0)
log_z_space = log10(z_space_magnitude + 1e-3);

% "Limit the display range by thresholding too large/too small values"
min_value = min(log_z_space(:));
max_value = max(log_z_space(:));
disp(min_value);
disp(max_value);
percent_cap = 0.05;
lower_threshold = min_value + percent_cap * (max_value - min_value); 
upper_threshold = max_value - percent_cap * (max_value - min_value);  
log_z_space_clipped = log_z_space;
log_z_space_clipped(log_z_space < lower_threshold) = lower_threshold;
log_z_space_clipped(log_z_space > upper_threshold) = upper_threshold;

% Display B-scan
figure;
imagesc(log_z_space_clipped);
colormap(gray);
axis image
title('B-Scan of Human Eye');
xlabel('x-pixel');
ylabel('y-pixel');