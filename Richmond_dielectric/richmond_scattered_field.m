function Esca = richmond_scattered_field_xy(xobs, yobs, Xm, Ym, eps_vec, E_unknown, k0, a_eq)

j = 1i;
N = numel(Xm);

% Precompute constants
J1 = besselj(1, k0 * a_eq);
H1 = besselh(1, 2, k0 * a_eq);

% Self-term integral
SelfInt = 1i / 2 * (pi * a_eq * k0 * H1 - 2*j);

Esca = zeros(size(xobs));

for p = 1:numel(xobs)
    xo = xobs(p);
    yo = yobs(p);

    Es = 0;

    for n = 1:N
        dx = xo - Xm(n);
        dy = yo - Ym(n);
        R  = sqrt(dx^2 + dy^2);

        if R < a_eq
            % Observation point inside cell n
            Int = SelfInt;
        else
            % Mutual-cell approximation
            H0 = besselh(0, 2, k0 * R);
            Int = (1i * pi * a_eq * k0 / 2) * J1 * H0;
        end

        Es = Es + (eps_vec(n) - 1) * E_unknown(n) * Int;
    end

    Esca(p) = Es;
end
end
