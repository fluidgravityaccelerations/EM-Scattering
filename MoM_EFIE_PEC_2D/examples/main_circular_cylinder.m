% --- Method of moments - One point quadrature - Scattering by a circular
% cylinder

close all 
clear all
clc

epsilon0    = 8.85418781 * 1e-12;                                           % --- Free-space dielectric permittivity
mu0         = 4 * pi * 1e-7;                                                % --- Free-space magnetic permeability
zeta0       = sqrt(mu0 / epsilon0);                                         % --- Free-space wave impendance 
c           = 1 / sqrt(epsilon0 * mu0);                                     % --- Wavespeed in vacuum

f           = 3e9;                                                          % --- Operating frequency
omega       = 2 * pi * f;                                                   % --- Operating angular frequency
beta        = omega / c;                                                    % --- Wavenumber
lambda      = 2 * pi / beta;                                                % --- Wavelength

a           = 5 * lambda;                                                             % --- Cylinder radius  

NN          = ceil(2 * pi * a / (lambda / 10));                             % --- Number of discretization intervals of the cylinder's contour
                                            
afar        = 10 * a * a / lambda;                                          % --- Far-field calculation distance

E0 = 1;                                                                     % --- Amplitude of the impinging wave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION OF THE CYLINDER SURFACE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x, y, xp, yp, tp, Deltal] = discretize_cylinder(a, NN);

%%%%%%%%%%%%%%%%%%%%%%%
% COEFFICIENTS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%
A = build_mom_matrix(xp, yp, beta, omega, mu0, Deltal);

%%%%%%%%%%%%%%%%%%%
% IMPINGING FIELD %
%%%%%%%%%%%%%%%%%%%
Ei = build_incident_field(xp, beta, E0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUTION - CURRENT RECONSTRUCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = linsolve(A, Ei);

Jref = compute_reference_solution(tp, beta, a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCATTERED FAR-FIELD CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Esfar = compute_far_field(tp, J, beta, omega, mu0, a, Deltal, afar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION OF THE EXTERIOR REGION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ae          = 3 * a;                                                        % --- Semi-extent of "exterior" region along x-axis
be          = 3 * a;                                                        % --- Semi-extent of "exterior" region along y-axis
Nx          = ceil(ae / (lambda / 8));                                      % --- Number of discretization points of the "exterior" region along x
Ny          = ceil(ae / (lambda / 8));                                      % --- Number of discretization points of the "exterior" region along y
[Xe, Ye]    = meshgrid(linspace(-ae, ae, Nx), linspace(-be, be, Ny));
indinterior = find(sqrt(Xe.^2 + Ye.^2) <= a);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEAR-FIELD CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%
[Es, E] = compute_near_field(Xe, Ye, xp, yp, J, beta, omega, mu0, Deltal, E0, a);

%%%%%%%%%
% PLOTS %
%%%%%%%%%
figure(1);
plot(tp, abs(J), 'k');
title('Amplitude of the induced current density')

figure(2);
plot(tp, abs(Esfar), 'k');
title('Amplitude of the far-field')

figure(3)
surf(Xe, Ye, abs(Es)); 
view(2), axis([-ae ae -be be]), axis equal, shading interp, colorbar
title('Amplitude of the scattered field')

figure(4)
surf(Xe, Ye, abs(E)); 
view(2), axis([-ae ae -be be]), axis equal, shading interp, colorbar
title('Amplitude of the total field')

%%%%%%%%%%%%%
% ANIMATION %
%%%%%%%%%%%%%
animate_field(Xe, Ye, E, f, ae, be);



