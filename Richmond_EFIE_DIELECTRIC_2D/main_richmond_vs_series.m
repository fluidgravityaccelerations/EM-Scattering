clear; clc; close all;

%% Physical parameters
f0      = 3e9;
c0      = 3e8;
lambda0 = c0 / f0;
k0      = 2*pi / lambda0;
eps_r   = 2.5;
theta_i = 0;

%% Geometry
a_cyl = 1.0 * lambda0;   % cylinder radius

%% Richmond grid
Lx   = 4*lambda0;
Ly   = 4*lambda0;
Nx   = ceil(Lx / (lambda0 / 30));
Ny   = ceil(Ly / (lambda0 / 30));

[Xc, Yc, dx, dy] = richmond_make_grid(Lx, Ly, Nx, Ny);
Acell = dx * dy;
a_eq  = sqrt(Acell / pi);

eps_map = ones(size(Xc));
inside  = (Xc.^2 + Yc.^2) <= a_cyl^2;
eps_map(inside) = eps_r;

mask        = (eps_map ~= 1);
idx         = find(mask);
Xm          = Xc(idx);
Ym          = Yc(idx);
eps_vec     = eps_map(idx);
Nunk        = numel(idx);
fprintf('Richmond unknowns: %d\n', Nunk);

%% Assemble and solve Richmond system
C           = richmond_assemble_matrix(Xm, Ym, eps_vec, k0, a_eq);
C           = C + eye(size(C));
Einc_full   = richmond_plane_wave_incident(Xc, Yc, k0, theta_i);
b           = Einc_full(idx);
E_unknown   = C \ b;
%% Reconstruct total field on full grid
E_tot       = Einc_full;     % start with incident everywhere
E_tot(idx)  = E_unknown;  % replace inside dielectric

%% Scattered field from Richmond on observation circle
Nobs   = 360;
phi_obs = linspace(0, 2*pi, Nobs);
Robs    = 3.0 * lambda0;

xobs = Robs * cos(phi_obs);
yobs = Robs * sin(phi_obs);

Esca_rich = richmond_scattered_field(xobs, yobs, Xm, Ym, eps_vec, E_unknown, k0, a_eq);

%% Cylindrical-harmonics reference
Nmax = 2 * ceil(k0 * a_cyl);   % truncation; increase if needed
Esca_series = tmz_dielectric_cylinder_series(phi_obs, Robs, k0, eps_r, a_cyl, theta_i, Nmax);

err_rel = abs(Esca_rich - Esca_series) ./ max(abs(Esca_series));
fprintf('Max relative error (magnitude, normalized): %.2e\n', max(err_rel));

%% Plots
figure;
plot(phi_obs*180/pi, 20*log10(abs(Esca_series)), 'k-', 'LineWidth', 1.5); hold on;
plot(phi_obs*180/pi, 20*log10(abs(Esca_rich)), 'r--', 'LineWidth', 1.2);
xlabel('\phi [deg]'); ylabel('|E_s| (dB, normalized)');
legend('Series (reference)', 'Richmond VIE', 'Location', 'Best');
title('Scattered field pattern: Richmond vs. cylindrical series');
grid on;

figure;
plot(phi_obs*180/pi, angle(Esca_series) * 180 / pi, 'k-', 'LineWidth', 1.5); hold on;
plot(phi_obs*180/pi, angle(Esca_rich) * 180 / pi, 'r--', 'LineWidth', 1.2);
xlabel('\phi [deg]'); ylabel('|E_s| (dB, normalized)');
legend('Series (reference)', 'Richmond VIE', 'Location', 'Best');
title('Scattered field pattern: Richmond vs. cylindrical series');
grid on;

figure;
plot(phi_obs*180/pi, err_rel, 'b-','LineWidth',1.2);
xlabel('\phi [deg]'); ylabel('Relative error');
title('Pointwise relative error |E_{Rich}/E_{series} - 1|');
grid on;

