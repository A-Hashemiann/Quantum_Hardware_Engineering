% Lab on Band Diagram: 1D Periodic Crystal
% Parameters
a = 1;      % Lattice constant
V0 = 1;     % Potential strength
N = 20;     % Basis: m = -N to N
m = -N:N;
G = (2*pi/a) * m';  % Column vector
num_basis = length(G);

% k-points in Brillouin zone
nk = 100;
klist = linspace(-pi/a, pi/a, nk);
E = zeros(4, nk);  % First 4 bands

% Compute bands
for ik = 1:nk
    k = klist(ik);
    q = k + G;
    H = diag(q.^2 / 2);
    H = H + diag(V0 * ones(num_basis-1,1), 1);
    H = H + diag(V0 * ones(num_basis-1,1), -1);
    [Vec, Val] = eig(H);
    Val = diag(Val);
    [Val, idx] = sort(Val);
    E(:, ik) = Val(1:4);
end

% Plot band diagram
figure(1);
plot(klist / (pi/a), E', 'LineWidth', 1.5);
xlabel('k (\pi/a)');
ylabel('Energy (atomic units)');
title('First Four Energy Bands');
grid on;
legend('Band 1', 'Band 2', 'Band 3', 'Band 4');

% Wavefunctions for specific k and bands
x = linspace(0, a, 100);
k_values = [0, pi/a];  % Small and large k
band_indices = [1, 2];

for ik = 1:2
    k = k_values(ik);
    q = k + G;
    H = diag(q.^2 / 2);
    H = H + diag(V0 * ones(num_basis-1,1), 1);
    H = H + diag(V0 * ones(num_basis-1,1), -1);
    [Vec, Val] = eig(H);
    Val = diag(Val);
    [~, idx] = sort(Val);
    Vec = Vec(:, idx);
    
    for ib = 1:2
        band = band_indices(ib);
        c = Vec(:, band);
        u = zeros(size(x));
        for im = 1:num_basis
            u = u + c(im) * exp(1i * G(im) * x);
        end
        % Normalize: int_0^a |u|^2 dx ≈ a (approximate with trapz)
        norm_integral = trapz(x, abs(u).^2);
        u = u * sqrt(a / (norm_integral + eps));  % Avoid div by zero
        
        % Plot real part (adjust phase if needed to make real)
        figure(1 + ik*2 + ib);
        plot(x, real(u), 'LineWidth', 1.5);
        xlabel('x');
        ylabel('Re(u_k(x))');
        title(sprintf('Wavefunction: k=%.2f \\pi/a, Band %d', k/(pi/a), band));
        grid on;
    end
end