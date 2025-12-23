% constant material properties
epsilon_0 = 8.85e-12;
mu_0 = 4e-7 * pi;

% indices of refraction - fixed values
n_i = 1;
n_s = 1.5;
n_l = 1.48;
n_h = 2.4;

% layer lengths according to 4.1 q1
d_h = 92.2e-9;
d_l = 163.2e-9;

% permittivities - fixed values
epsilon_i = n_i ^ 2 * epsilon_0;
epsilon_s = n_s ^ 2 * epsilon_0;
epsilon_l = n_l ^ 2 * epsilon_0;
epsilon_h = n_h ^ 2 * epsilon_0;

% incident angles [rad] - fixed values
theta_i = pi/4;
theta_s = 0.49;
theta_l = 0.5;
theta_h = 0.3;

% range of wavelengths
samples = 700;
lambda = 1e-9 * linspace(500, 1200, samples);

% Impedances - TM
Z_i_TM = sqrt(mu_0 / epsilon_i) * (cos(theta_i));
Z_s_TM = sqrt(mu_0 / epsilon_s) * (cos(theta_s));
Z_l_TM = sqrt(mu_0 / epsilon_l) * (cos(theta_l));
Z_h_TM = sqrt(mu_0 / epsilon_h) * (cos(theta_h));

% adjustable number of layers
N = 20;

% Initialize the reflection array
reflection = zeros(1, samples);

function matrix = propagation_matrix(Z_l_plus_1, Z_l, kz_l_plus_1, d_l_plus_1, d_l)
    Gamma_l_1_l = (Z_l - Z_l_plus_1) / (Z_l + Z_l_plus_1);
    P = Z_l_plus_1 / Z_l;
    M11 = exp(-1j * kz_l_plus_1 * (d_l_plus_1 - d_l));
    M12 = Gamma_l_1_l * exp(-1j * kz_l_plus_1 * (d_l_plus_1 - d_l));
    M21 = Gamma_l_1_l * exp(+1j * kz_l_plus_1 * (d_l_plus_1 - d_l));
    M22 = exp(+1j * kz_l_plus_1 * (d_l_plus_1 - d_l));
    constants = (1 / 2) * (1 + P);
    matrix = constants * [M11, M12; M21, M22];
end


% Loop through each wavelength
for i = 1:samples
    ki = 2 * pi / lambda(i);
    kz_h = 2 * pi / (lambda(i) / n_h) * cos(theta_h);
    kz_l = 2 * pi / (lambda(i) / n_l) * cos(theta_l);
    kz_s = 2 * pi / (lambda(i) / n_s) * cos(theta_s);
    
    V_hi = propagation_matrix(Z_h_TM, Z_i_TM, kz_h, d_h, 0);
    V_lh = propagation_matrix(Z_l_TM, Z_h_TM, kz_l, d_l, 0);
    V_hl = propagation_matrix(Z_h_TM, Z_l_TM, kz_h, d_h, 0);
    V_sh = propagation_matrix(Z_s_TM, Z_h_TM, kz_s, 0, d_h);
    
    full_matrix = V_sh * (V_hl * V_lh) ^ N * V_hi;
    reflection(i) = -full_matrix(2, 1) / full_matrix(2, 2);
end

% Plot Power Reflected vs. Wavelength
figure;
subplot(2, 1, 1);
plot(lambda * 1e9, abs(reflection).^2); % Convert lambda to nm for plotting
title('Power Reflected vs. Wavelength (TM) (N = 20)');
xlabel('Wavelength (nm)');
ylabel('Power Reflected');
grid on;

subplot(2, 1, 2);
plot(lambda * 1e9, angle(reflection)); % Convert lambda to nm for plotting
title('Reflection Phase vs. Wavelength (TM) (N = 20)');
xlabel('Wavelength (nm)');
ylabel('Reflection Phase (radians)');
grid on;