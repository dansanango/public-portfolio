% constant material properties
epsilon_0 = 8.85e-12;
mu_0 = 4e-7 * pi;

% indices of refraction-fixed values
n_i = 1;
n_s = 1.5;
n_l = 1.48;
n_h = 2.4;

% permittivities-fixed values
epsilon_i = n_i ^ 2 * epsilon_0;
epsilon_s = n_s ^ 2 * epsilon_0;
epsilon_l = n_l ^ 2 * epsilon_0;
epsilon_h = n_h ^ 2 * epsilon_0;

% incident angles [rad]-fixed values
theta_i = pi/4;
theta_s = 0.49;
theta_l = 0.5;
theta_h = 0.3;

% Impedances
Z_i_TE = sqrt(mu_0 / epsilon_i) * (1 / cos(theta_i));
Z_s_TE = sqrt(mu_0 / epsilon_s) * (1 / cos(theta_s));
Z_l_TE = sqrt(mu_0 / epsilon_l) * (1 / cos(theta_l));
Z_h_TE = sqrt(mu_0 / epsilon_h) * (1 / cos(theta_h));

% Define the Z_z function
function Z_z = Z_z(N, Z_s_TE, Z_l_TE, Z_h_TE)
    Z_z = (1/Z_s_TE) * (Z_h_TE/Z_l_TE)^(2*N) * (Z_h_TE)^2;
end

% Define the Gamma function
function Gamma = Gamma(N, Z_s_TE, Z_l_TE, Z_h_TE, Z_i_TE)
    Z_z_val = Z_z(N, Z_s_TE, Z_l_TE, Z_h_TE);
    Gamma = (Z_z_val - Z_i_TE) / (Z_z_val + Z_i_TE);
end

% Initialize N and calculate the power reflection
N = 0;
Power_R = Gamma(N, Z_s_TE, Z_l_TE, Z_h_TE, Z_i_TE)^2;

% Loop until Power_R exceeds 99.9% = 0.999
while Power_R <= 0.999
    N = N + 1;
    Power_R = Gamma(N, Z_s_TE, Z_l_TE, Z_h_TE, Z_i_TE)^2;
end

% Display the final values
disp("N:");
disp(N);
disp("Power Reflected:");
disp(Gamma(N, Z_s_TE, Z_l_TE, Z_h_TE, Z_i_TE)^2);
disp("Reflection Coefficient (should be about -1):");
disp(Gamma(N, Z_s_TE, Z_l_TE, Z_h_TE, Z_i_TE));
