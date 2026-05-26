function A = build_mom_matrix(xp, yp, beta, omega, mu0, Deltal)
%BUILD_MOM_MATRIX  Constructs the MoM coefficient matrix for a 2D PEC cylinder.
%
%   A = build_mom_matrix(xp, yp, beta, omega, mu0, Deltal)
%
%   Inputs:
%       xp, yp   - Matching point coordinates (size NN)
%       beta     - Wavenumber (omega / c)
%       omega    - Angular frequency
%       mu0      - Free-space magnetic permeability
%       Deltal   - Effective segment length
%
%   Output:
%       A        - NN x NN Method of Moments coefficient matrix
%
%   Notes:
%   - The kernel is the 2D Green's function:
%
%         G(R) = (omega * mu0 / 4) * H0^(2)(beta * R)
%
%   - The diagonal term is replaced by the analytic limit of the
%     Hankel function singularity:
%
%         H0^(2)(kR) ~ (2j/pi) * [ ln(kR/2) + gamma ]
%
%     where gamma is the Euler–Mascheroni constant.
%
%   - In the integrated form over a segment of length Δl, the singular
%     term becomes:
%
%         log(e^γ * k * Δl / 4)
%
%     with e^γ ≈ 1.781072.
%

    NN = length(xp);

    % Pairwise distances between matching points
    Xp = xp(:) * ones(1, NN) - ones(NN, 1) * xp(:).';
    Yp = yp(:) * ones(1, NN) - ones(NN, 1) * yp(:).';
    R  = sqrt(Xp.^2 + Yp.^2);

    % Green's function (non-diagonal part)
    A = (omega * mu0 / 4) * besselh(0, 2, beta * R) * Deltal;

    % Euler–Mascheroni constant in exponential form
    % e^γ = 1.781072
    gamma_exp = 1.781072;

    % Analytic diagonal term
    diag_term = (omega * mu0 / 4) * Deltal * ...
                (-2i / pi) * ( log(gamma_exp * beta * Deltal / 4) + 1i*pi/2 - 1 );

    % Replace diagonal entries
    A(1:NN+1:end) = diag_term;

end
