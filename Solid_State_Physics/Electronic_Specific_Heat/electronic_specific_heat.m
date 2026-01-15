clc; clear; close all;


%% Physical constants (SI units)
hbar = 1.054e-34; % Planck constant / 2pi [J s]
m = 9.11e-31; % electron mass [kg]
kB = 1.38e-23; % Boltzmann constant [J/K]


%% Parameters
Ef = 10 * 1.602e-19; % Fermi energy [J]
T = linspace(1,300,200);


%% Energy grid
E = linspace(0, 2*Ef, 5000);
dE = E(2) - E(1);


%% Density of states
D = (1/(2*pi^2)) * (2*m/hbar^2)^(3/2) * sqrt(E);


U = zeros(size(T));


%% Compute internal energy
for i = 1:length(T)
f = 1 ./ (exp((E - Ef) ./ (kB*T(i))) + 1);
integrand = D .* E .* f;
U(i) = trapz(E, integrand);
end


%% Plot internal energy vs temperature
figure;
plot(T, U, 'LineWidth', 2);
xlabel('Temperature (K)');
ylabel('Internal Energy (J)');
title('Electronic Internal Energy vs Temperature');
grid on;


%% Electronic specific heat (numerical derivative)
ce = gradient(U, T);


%% Low-temperature linear fit
idx_lowT = T < 50;
coeffs = polyfit(T(idx_lowT), ce(idx_lowT), 1);


figure;
plot(T, ce, 'LineWidth', 2);
hold on;
plot(T(idx_lowT), polyval(coeffs, T(idx_lowT)), '--r', 'LineWidth', 2);
xlabel('Temperature (K)');
ylabel('c_e (J/K)');
title('Electronic Specific Heat');
legend('Numerical', 'Low-T linear fit');
grid on;