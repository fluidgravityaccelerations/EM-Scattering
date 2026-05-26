function J_exact = current_circle(phi, a, k, eta0)
%current_circle Exact TMz surface current on a circular PEC cylinder.
%
%   J_exact = current_circle(phi, a, k, eta0)
%
%   phi  : azimuthal angle (vector)
%   a    : radius
%   k    : wavenumber
%   eta0 : intrinsic impedance
%
%   Plane wave: E_inc = E0 * exp(-j*k*x)
%   Time convention: exp(-j*omega*t), outgoing waves ~ H_n^{(2)}.

    E0   = 1;
    Nmax = ceil(k*a + 10);
    J_exact = zeros(size(phi));

    % Series: sum_{n=-∞}^{∞} j^{-n} e^{j n phi} / H_n^{(2)}(k a)
    % Implemented as i^(-n) e^{i n phi} / H_n^{(2)}(k a)
    for n = -Nmax:Nmax
        J_exact = J_exact + (1i^(-n)) * exp(1i*n*phi) / besselh(n,2,k*a);
    end

    % Prefactor: 2 E0 / (pi a omega mu0) = 2 E0 / (pi a k eta0)
    J_exact = (2*E0 / (pi * a * k * eta0)) * J_exact;
end
