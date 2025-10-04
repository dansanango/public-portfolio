% Sellmeier coefficients
A_o = 2.7359; B_o = 0.01878; C_o = 0.01822; D_o = 0.01354;
A_e = 2.3753; B_e = 0.01224; C_e = 0.01667; D_e = 0.01516;

% Threshold 
threshold = 1e-4;

% Define wavelength
lambda_range = linspace(0.2, 2.5, 1000); 
% Define angles
theta_range = linspace(0, pi/2, 1000);

% Sellmeier for no and ne
no_values = sqrt(A_o + B_o ./ (lambda_range.^2 - C_o) - D_o * lambda_range.^2);
ne_values = sqrt(A_e + B_e ./ (lambda_range.^2 - C_e) - D_e * lambda_range.^2);

% Lists to store valid (lambda, theta) pairs
lambda_valid = [];
theta_valid = [];
ne_theta_values = [];
ne_valid = [];
no_valid = [];

% Iterate over all phase matching angles
for i = 1:length(theta_range)
    theta = theta_range(i);
    
    % Iterate over each lambda in the discrete list
    for j = 1:length(lambda_range)
        lambda = lambda_range(j);
        lambda_half = lambda / 2;

        % Ensure lambda_half is within the valid range before interpolation
        if lambda_half >= min(lambda_range)
            no_lambda = no_values(j);
            ne_lambda = ne_values(j);
            
            % Interpolate n_e at lambda/2
            ne_lambda_half = interp1(lambda_range, ne_values, lambda_half, 'linear', 'extrap');
            no_lambda_half = interp1(lambda_range, no_values, lambda_half, 'linear', 'extrap');
            % Compute extraordinary refractive index at theta
            ne_theta = 1 / sqrt((cos(theta) / no_lambda_half)^2 + (sin(theta) / ne_lambda_half)^2);
            
            % Check phase matching condition
            if abs(no_lambda - ne_theta) < threshold
                % a bunch of ugly code to get all good data in a list
                lambda_valid = [lambda_valid, lambda]; %#ok<AGROW>
                theta_valid = [theta_valid, theta]; %#ok<AGROW>
                ne_theta_values = [ne_theta_values, ne_theta]; %#ok<AGROW>
                ne_valid = [ne_valid, ne_lambda];%#ok<AGROW>
                no_valid = [no_valid, no_lambda];%#ok<AGROW>
                
            end
        end
    end
end

theta_valid_deg = rad2deg(theta_valid);

figure;
plot(theta_valid_deg, lambda_valid, 'b.', 'MarkerSize', 10);
xlabel('Phase Matching Angle (degrees)');
ylabel('Phase Matching Wavelength (μm)');
title('Phase Matching Wavelength vs. Angle (Type I)');
grid on;