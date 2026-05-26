function Jref = compute_reference_solution(tp, beta, a)
%COMPUTE_REFERENCE_SOLUTION  Computes the analytical current on a PEC cylinder.
%
%   Jref = compute_reference_solution(tp, beta, a)
%
%   Inputs:
%       tp    - Angular positions of the matching points (size NN)
%       beta  - Wavenumber
%       a     - Cylinder radius
%
%   Output:
%       Jref  - Analytical surface current density (size NN)
%
%   Notes:
%   - This function computes the exact solution for the induced current
%     on a perfectly conducting circular cylinder illuminated by a
%     plane wave. The analytical expression is:
%
%         J(θ) = Σ_{n=-∞}^{∞} j^{-n} * exp(j n θ) / H_n^{(2)}(β a)
%
%   - The summation is truncated at:
%
%         n_max = ceil(2 * β * a)
%
%     which is a standard rule ensuring convergence.
%
%   - This reference solution is useful for validating the MoM numerical
%     solution and for error analysis.
%

    % Truncation index
    nmax = ceil(2 * beta * a);

    % Initialize reference current
    Jref = zeros(size(tp));

    % Summation over cylindrical harmonics
    for n = -nmax : nmax
        Jref = Jref + 1i^(-n) * exp(1i * n * tp) ./ besselh(n, 2, beta * a);
    end

end
