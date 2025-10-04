% Constants
B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;
C1 = 0.0684043^2;
C2 = 0.1162414^2;
C3 = 9.896161^2;
samples = 1e3;
c = 3e8;
um_to_m = 1e-6;
m_to_um = 1e6;

% Define lambda in [m] and [um]
lambda_m = linspace(0.6, 1.6, samples) * um_to_m;
lambda_um = lambda_m * m_to_um;

% index of refraction
n = sqrt(1 + (B1 .* lambda_um.^2 ./ (lambda_um.^2 - C1)) + (B2 .* lambda_um.^2 ./ (lambda_um.^2 - C2)) + (B3 .* lambda_um.^2 ./ (lambda_um.^2 - C3)));

% dn/dlambda
dn_dlambda = gradient(n, lambda_m);

% dbeta/domega
dbeta_domega = (1/c)*(n - lambda_m.*dn_dlambda);

% dbeta/domega/dlambda (dispersion)
d_dbeta_domega_dlambda = gradient(dbeta_domega, lambda_m);

% Convert to dispersion in ps/nm/km
dispersion_ps_nm_km = d_dbeta_domega_dlambda * 1e6;

% Plot
figure;
plot(lambda_um, dispersion_ps_nm_km, 'r', 'LineWidth', 2);
xlabel('Wavelength (microns)');
ylabel('Dispersion (ps/nm/km)');
title('Dispersion vs. Wavelength');
grid on;
