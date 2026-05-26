function Esca = tmz_dielectric_cylinder_series(phi_obs, Robs, k0, eps_r, a_cyl, theta_i, Nmax)
%TMZ_DIELECTRIC_CYLINDER_SERIES Cylindrical-harmonics solution (TMz).
%
%   Esca = tmz_dielectric_cylinder_series(phi_obs, Robs, k0, eps_r, a_cyl, theta_i, Nmax)
%
%   phi_obs : observation angles (vector)
%   Robs    : observation radius (scalar, > a_cyl)
%   k0      : free-space wavenumber
%   eps_r   : relative permittivity of cylinder
%   a_cyl   : cylinder radius
%   theta_i : incidence angle of plane wave (from +x axis)
%   Nmax    : truncation index for cylindrical harmonics
%
%   Esca    : scattered Ez at (Robs, phi_obs)

j   = 1i;
k1  = k0 * sqrt(eps_r);   % wavenumber inside cylinder

phi_obs = phi_obs(:).';   % row vector
Esca    = zeros(size(phi_obs));

for n = -Nmax:Nmax
    % Incident field expansion coefficient (TMz, plane wave)
    % E_inc = sum_n j^(-n) J_n(k0 r) e^{j n (phi - theta_i)}
    A_n = 1i^(-n) * exp(-1i * n * theta_i);

    % Boundary conditions at r = a_cyl:
    % continuity of Ez and (1/mu) dEz/dr (mu same -> continuity of dEz/dr)
    Jn_k0a  = besselj(n, k0 * a_cyl);
    Jn_k1a  = besselj(n, k1 * a_cyl);
    Jn_p_k0a = besseljp(n, k0 * a_cyl);   % derivative wrt argument
    Jn_p_k1a = besseljp(n, k1 * a_cyl);
    Hn_k0a  = besselh(n, 2, k0 * a_cyl);
    Hn_p_k0a = besselhp(n, 2, k0 * a_cyl);

    % Solve for scattered coefficient B_n (outside)
    % A_n J_n(k0 a) + B_n H_n(k0 a) = C_n J_n(k1 a)
    % k0 [A_n J_n'(k0 a) + B_n H_n'(k0 a)] = k1 C_n J_n'(k1 a)
    %
    % Eliminate C_n, solve for B_n:
    num = A_n * ( Jn_k0a * k1 * Jn_p_k1a - k0 * Jn_p_k0a * Jn_k1a );
    % den = Hn_k0a * k0 * Jn_p_k1a - Hn_p_k0a * k1 * Jn_k1a;
    den = Hn_k0a * k1 * Jn_p_k1a - Hn_p_k0a * k0 * Jn_k1a;
    B_n = num / den;

    % Contribution to scattered field at (Robs, phi_obs)
    Hn_k0R = besselh(n, 2, k0 * Robs);
    Esca = Esca + B_n * Hn_k0R .* exp(1i * n * phi_obs);
end
end

% --- helper: derivative of Bessel J_n ---
function Jp = besseljp(n, z)
% derivative wrt argument: J_n'(z) = 0.5*(J_{n-1}(z) - J_{n+1}(z))
Jp = 0.5 * (besselj(n-1, z) - besselj(n+1, z));
end

% --- helper: derivative of Hankel H_n^{(2)} ---
function Hp = besselhp(n, kind, z)
% derivative wrt argument: H_n'(z) = 0.5*(H_{n-1}(z) - H_{n+1}(z))
Hn_1 = besselh(n-1, kind, z);
Hn1  = besselh(n+1, kind, z);
Hp   = 0.5 * (Hn_1 - Hn1);
end
