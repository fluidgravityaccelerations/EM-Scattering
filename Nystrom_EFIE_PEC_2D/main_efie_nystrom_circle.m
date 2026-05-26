% -------------------------------------------------------------------------
% main_efie_nystrom_circle.m
%
% EFIE Nyström solver for 2D scattering by a circular/elliptical PEC cylinder
% Colton–Kress trigonometric quadrature + Martensen–Kussmaul splitting
% EFIE (TM polarization): only the single-layer operator (M1, M2)
%
% Time convention: exp(j*omega*t)
% Green's function: G = (j/4) H_0^{(2)}(kR)
% -------------------------------------------------------------------------

close all; clear all; clc;

% ------------------------------
% Electromagnetic parameters
% ------------------------------
epsilon0 = 8.8541878128e-12;
mu0      = 4*pi*1e-7;
c        = 1/sqrt(epsilon0*mu0);
lambda   = 1;
k        = 2*pi/lambda;
omega    = k * c;
eta0     = sqrt(mu0/epsilon0);   % intrinsic impedance of free space

% ------------------------------
% Geometry (circle if a = b)
% ------------------------------
a = 2*lambda;   % semi-axis along x (radius if a = b)
b = 2*lambda;   % semi-axis along y

% ------------------------------
% Nyström discretization
% ------------------------------
N  = 150;                                   % number of intervals
tn = pi*(2*(0:2*N-1)+1)/(2*N);              % 2N midpoint nodes in [0, 2*pi)

% Parametrization of the boundary
Xn  = a*cos(tn);
Yn  = b*sin(tn);

% First derivatives (tangent vector)
dXn = -a*sin(tn);
dYn =  b*cos(tn);

% Arc-length metric |x'(t)|
ds = sqrt(dXn.^2 + dYn.^2);

% ------------------------------
% Allocate matrices
% ------------------------------
M1 = zeros(2*N);
M2 = zeros(2*N);
R  = zeros(2*N);
A  = zeros(2*N);

% Euler–Mascheroni constant
gamma = 0.5772156649015328606;

% ------------------------------
% Build EFIE operator (single-layer only)
% ------------------------------
for i = 1:2*N
    for j = 1:2*N
        
        Rij = sqrt((Xn(i)-Xn(j))^2 + (Yn(i)-Yn(j))^2);  % distance
        dRi = ds(i);
        dRj = ds(j);

        % Colton–Kress trigonometric weight for the log kernel
        R(i,j) = RColtonKress(N, tn(i), tn(j));

        % --- M1: coefficient of the logarithmic singularity
        % M1(t_i, t_j) = -(1/(2*pi)) J_0(k R_ij) |x'(t_j)|
        M1(i,j) = -(1/(2*pi)) * besselj(0, k*Rij) * dRj;

        % --- M2: regular part (Martensen–Kussmaul splitting)
        if i == j
            % Diagonal limit (self-term)
            % M2(t_i, t_i) = (-j/2 - gamma/pi - (1/pi) log((k/2)|x'(t_i)|)) |x'(t_i)|
            M2(i,j) = (1i/2 - gamma/pi - (1/pi)*log((k/2)*dRi)) * dRi;
        else
            % Off-diagonal regular part
            logterm = log(4*(sin((tn(i)-tn(j))/2))^2);
            % M2(t_i, t_j) = (-j/2) H_0^{(2)}(k R_ij) |x'(t_j)| - M1 log(...)
            M2(i,j) = (1i/2)*besselh(0,1,k*Rij)*dRj - M1(i,j)*logterm;
        end

        % --- EFIE operator: j*omega*mu0 * S[J]
        % Discrete Nyström: A_ij = j*omega*mu0 [ R_ij M1_ij + (pi/N) M2_ij ]
        % 0.5 factor is a convention-dependent scaling (kept from derivation)
        A(i,j) = 0.5 * 1i*omega*mu0 * conj( R(i,j)*M1(i,j) + (pi/N)*M2(i,j) );
    end
end

% ------------------------------
% Right-hand side (incident field on boundary)
% Plane wave: E_inc = exp(-j*k*x), propagating along +x
% ------------------------------
g = exp(-1i*k*Xn).';   % column vector

% EFIE on PEC: E_inc + j*omega*mu0 S[J] = 0  =>  A J = -g
J = A \ g;

% ------------------------------
% Reference current for circular cylinder (Balanis-like)
% ------------------------------
J_ref = current_circle(tn, a, k, eta0);

% ------------------------------
% Plot: magnitude comparison
% ------------------------------
figure(1);
plot(tn, abs(J), 'k', 'LineWidth', 1.2); hold on;
plot(tn, abs(J_ref), 'r--', 'LineWidth', 1.2);
xlabel('t (rad)');
ylabel('|J(t)|');
legend('Nyström', 'Reference (series)', 'Location', 'best');
title('Induced surface current density - Magnitude');
grid on;

% ------------------------------
% Plot: phase comparison
% ------------------------------
figure(2);
plot(tn, angle(J), 'k', 'LineWidth', 1.2); hold on;
plot(tn, angle(J_ref), 'r--', 'LineWidth', 1.2);
xlabel('t (rad)');
ylabel('Phase[J(t)] (rad)');
legend('Nyström', 'Reference (series)', 'Location', 'best');
title('Induced surface current density - Phase');
grid on;

% ------------------------------
% Far-field pattern from Nyström current
% ------------------------------
Nphi   = 360;
phi_ff = linspace(0, 2*pi, Nphi+1); 
phi_ff = phi_ff(1:end-1);          % avoid duplicate 2*pi

F = zeros(size(phi_ff));           % angular far-field pattern

for m = 1:numel(phi_ff)
    phi = phi_ff(m);
    phase = k*(Xn*cos(phi) + Yn*sin(phi));   % dot product r_hat · r_n
    F(m)  = 1i*omega*mu0/4 * sum( exp(1i*phase) .* J.' .* ds );
end

figure(3);
plot(phi_ff, abs(F), 'k', 'LineWidth', 1.2);
xlabel('\phi (rad)');
ylabel('|F(\phi)|');
title('Far-field pattern from Nyström current');
grid on;

% ------------------------------
% Reference far-field (cylindrical Mie series)
% ------------------------------
F_ref = farfield_circle(phi_ff, a, k);

figure(4);
plot(phi_ff, abs(F)/max(abs(F)), 'k', 'LineWidth', 1.2); hold on;
plot(phi_ff, abs(F_ref)/max(abs(F_ref)), 'r--', 'LineWidth', 1.2);
xlabel('\phi (rad)');
ylabel('Normalized |F(\phi)|');
legend('Nyström', 'Reference (series)', 'Location', 'best');
title('Far-field pattern comparison');
grid on;
