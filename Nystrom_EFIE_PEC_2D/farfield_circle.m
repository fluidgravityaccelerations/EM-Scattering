function Fref = farfield_circle(phi, a, k)
%farfield_circle Reference far-field (TMz) for a circular PEC cylinder.
%
%   Fref = farfield_circle(phi, a, k)
%
%   phi : observation angle (vector)
%   a   : radius
%   k   : wavenumber
%
%   Based on the cylindrical Mie-series:
%   E_z^s ~ -E0 sqrt(2j/(pi k)) e^{-jk rho}/sqrt(rho)
%           * sum_n [ J_n(ka)/H_n^{(2)}(ka) e^{j n phi} ].

    E0 = 1;
    Nmax = ceil(k*a + 20);
    Fref = zeros(size(phi));

    for n = -Nmax:Nmax
        Fref = Fref + (besselj(n, k*a) / besselh(n, 2, k*a)) .* exp(1i*n*phi);
    end

    Fref = -E0 * sqrt(2*1i/(pi*k)) * Fref;
end
