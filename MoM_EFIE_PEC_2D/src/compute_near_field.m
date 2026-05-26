function [Es, E] = compute_near_field(Xe, Ye, xp, yp, J, beta, omega, mu0, Deltal, E0, a)
%COMPUTE_NEAR_FIELD  Computes the 2D near-field scattered and total field.
%
%   [Es, E] = compute_near_field(Xe, Ye, xp, yp, J, beta, omega, mu0, Deltal, E0, a)
%
%   Inputs:
%       Xe, Ye   - Evaluation grid coordinates (size Nx x Ny)
%       xp, yp   - Matching point coordinates (size NN)
%       J        - Solved surface current density (size NN)
%       beta     - Wavenumber
%       omega    - Angular frequency
%       mu0      - Free-space magnetic permeability
%       Deltal   - Effective segment length
%       E0       - Amplitude of the incident plane wave
%       a        - Cylinder radius (used to zero the interior field)
%
%   Outputs:
%       Es       - Scattered field on the grid (Nx x Ny)
%       E        - Total field (incident + scattered)
%
%   Notes:
%   - The scattered field is computed using the 2D Green's function:
%
%         G(R) = (ω μ0 / 4) * H0^(2)(β R)
%
%   - The computation is fully vectorized for efficiency.
%   - The total field is:
%
%         E = Es + E_inc
%
%     where the incident field is a plane wave exp(-j β x).
%

    % Grid size
    [Nx, Ny] = size(Xe);

    % Number of source panels
    NN = length(xp);

    % Preallocate scattered field
    Es = zeros(Nx, Ny);

    % Vectorized computation of distances:
    % For each panel k, compute R_k(x,y) = sqrt((x - xp(k))^2 + (y - yp(k))^2)
    for k = 1:NN
        Rk = sqrt((Xe - xp(k)).^2 + (Ye - yp(k)).^2);
        Es = Es - (omega * mu0 / 4) * Deltal * besselh(0, 2, beta * Rk) * J(k);
    end

    % Zero field inside the cylinder
    inside = (Xe.^2 + Ye.^2) <= a^2;
    Es(inside) = 0;

    % Incident plane wave
    Einc = E0 * exp(-1i * beta * Xe);

    % Total field
    E = Es + Einc;

end
