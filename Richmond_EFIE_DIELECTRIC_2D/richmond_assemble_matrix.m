function C = richmond_assemble_matrix(Xm, Ym, eps_vec, k0, a_eq)
%RICHMOND_ASSEMBLE_MATRIX Assemble Richmond VIE matrix for TMz scattering.
%
%   C = richmond_assemble_matrix(Xm, Ym, eps_vec, k0, a_eq)
%
%   Xm, Ym   : column vectors of cell-center coordinates (unknown cells)
%   eps_vec  : relative permittivity in each unknown cell
%   k0       : free-space wavenumber
%   a_eq     : equivalent circular radius of each cell (sqrt(area/pi))
%
%   C        : N x N complex matrix (Richmond's system matrix)

N = numel(Xm);
C = zeros(N, N);

% Precompute constants
j = 1i;

for m = 1:N
    for n = 1:N
        if m == n
            % Self-term (Richmond: singular integral treated analytically)
            % C_mm = 1 + (eps_m - 1) * j/2 * [pi*k0*a*H1^(2)(k0*a) - 2j]
            epsn = eps_vec(n);
            H1   = besselh(1, 2, k0 * a_eq);
            C(m,m) = (epsn - 1) * (j / 2) * (pi * k0 * a_eq * H1 - 2*j);
        else
            % Mutual interaction term
            % C_mn = j*pi*k0*a/2 * (eps_n - 1) * J1(k0*a) * H0^(2)(k0*R_mn)
            dx = Xm(m) - Xm(n);
            dy = Ym(m) - Ym(n);
            Rmn = sqrt(dx^2 + dy^2);

            epsn = eps_vec(n);
            J1   = besselj(1, k0 * a_eq);
            H0   = besselh(0, 2, k0 * Rmn);

            C(m,n) = (j * pi * k0 * a_eq / 2) * (epsn - 1) * J1 * H0;
        end
    end
end
end
