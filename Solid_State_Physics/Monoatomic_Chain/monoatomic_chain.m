clc; clear; close all;


N = 200; % Number of atoms
a = 1; % Lattice constant
M = 1; % Atomic mass
f = 1; % Force constant


omega_vals = linspace(0.01, 2*sqrt(f/M), 200);
k_vals = zeros(size(omega_vals));


for n = 1:length(omega_vals)
omega = omega_vals(n);
A = zeros(N);
for i = 1:N
A(i,i) = 2*f - M*omega^2;
if i > 1
A(i,i-1) = -f;
end
if i < N
A(i,i+1) = -f;
end
end
% Born–von Karman boundary conditions
A(1,N) = -f;
A(N,1) = -f;
[V, D] = eig(A);
[~, idx] = min(abs(diag(D)));
u = V(:, idx);
U = fft(u);
[~, idx_max] = max(abs(U(2:N/2)));
k_vals(n) = 2*pi*idx_max/(N*a);
end


% Adimensional plot
figure;
plot(k_vals*a, omega_vals/(2*sqrt(f/M)), 'b.', 'LineWidth', 1.5);
hold on;


ka = linspace(0, pi, 400);
omega_analytical = 2*sqrt(f/M)*abs(sin(ka/2));
plot(ka, omega_analytical/(2*sqrt(f/M)), 'r', 'LineWidth', 2);


xlabel('ka');
ylabel('\omega / (2\sqrt{f/M})');
legend('Numerical', 'Analytical');
grid on;

figure;
subplot(2,1,1);
plot(u, '-o');
title('Low ka mode');


grid on;


subplot(2,1,2);
plot(flipud(u), '-o');
title('High ka mode');


grid on;