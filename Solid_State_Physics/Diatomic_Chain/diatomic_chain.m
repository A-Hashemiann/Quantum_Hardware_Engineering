
%% -------------------------------
% Physical parameters
% -------------------------------
N = 1000;        % total number of atoms (must be even)
a = 1;           % lattice constant
f = 1;           % spring constant
M = 1;           % mass of atom 1
m = 2;           % mass of atom 2

%% -------------------------------
% Build mass matrix
% -------------------------------
Mass = zeros(N,1);
for i = 1:N
    if mod(i,2) == 1
        Mass(i) = M;
    else
        Mass(i) = m;
    end
end

%% -------------------------------
% Build force-constant matrix K
% -------------------------------
K = zeros(N);

for i = 1:N
    K(i,i) = 2*f;

    if i > 1
        K(i,i-1) = -f;
    end
    if i < N
        K(i,i+1) = -f;
    end
end

% Born–von Karman boundary conditions
K(1,N) = -f;
K(N,1) = -f;

%% -------------------------------
% Build dynamical matrix D = M^{-1} K
% -------------------------------
D = zeros(N);
for i = 1:N
    for j = 1:N
        D(i,j) = K(i,j) / sqrt(Mass(i)*Mass(j));
    end
end

%% -------------------------------
% Solve eigenvalue problem
% -------------------------------
[V, Lambda] = eig(D);

omega2 = diag(Lambda);
omega = sqrt(abs(omega2));   % angular frequencies

%% -------------------------------
% Extract wavevectors using FFT
% -------------------------------
k_vals = zeros(N,1);

for n = 1:N
    u = V(:,n);
    U = fft(u);
    [~, idx] = max(abs(U(2:N/2)));
    k_vals(n) = 2*pi*idx/(N*a);
end

%% -------------------------------
% Sort by k
% -------------------------------
[k_sorted, idx] = sort(k_vals);
omega_sorted = omega(idx);

%% -------------------------------
% Plot dispersion relation
% -------------------------------
figure;
scatter(k_sorted, omega_sorted, 10, 'filled');
xlabel('k');
ylabel('\omega');
title('Phonon dispersion relation – Diatomic chain');
grid on;
