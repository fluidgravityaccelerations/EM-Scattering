function Esfar = compute_far_field(tp, J, beta, omega, mu0, a, Deltal, afar)
%COMPUTE_FAR_FIELD  Computes the 2D far-field scattered by a PEC cylinder.
%
%   Esfar = compute_far_field(tp, J, beta, omega, mu0, a, Deltal, afar)
%
%   Inputs:
%       tp      - Angular positions of the matching points (size NN)
%       J       - Solved surface current density (size NN)
%       beta    - Wavenumber
%       omega   - Angular frequency
%       mu0     - Free-space magnetic permeability
%       a       - Cylinder radius
%       Deltal  - Effective segment length
%       afar    - Far-field observation distance
%
%   Output:
%       Esfar   - Far-field scattered electric field (size NN x 1)
%
%   Notes:
%   - The far-field expression for 2D scattering is:
%
%         E_s(θ) = C * Σ_n J_n * exp(j * β * a * cos(θ - θ_n))
%
%     where:
%         C = -(ω μ0 / 4) * sqrt(2 / (π r)) * Δl * exp(-j β r)
%
%   - The computation is fully vectorized using a kernel matrix K.
%

    NN = length(tp);

    % Constant prefactor for the far-field
    Cfar = - (omega * mu0 / 4) * sqrt(2 / (pi * afar)) * ...
            Deltal * exp(-1i * beta * afar);

    % Angular kernel matrix:
    % K(m,n) = exp( j * beta * a * cos( tp(n) - tp(m) ) )
    [TP_obs, TP_src] = meshgrid(tp, tp);
    K = exp(1i * beta * a .* cos(TP_src - TP_obs));   % NN x NN

    % Far-field as matrix–vector product
    Esfar = Cfar * (K * J);

end
