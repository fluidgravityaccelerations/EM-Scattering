function Einc = richmond_plane_wave_incident(X, Y, k0, theta_i)
%RICHMOND_PLANE_WAVE_INCIDENT TMz plane wave Ez^inc.
%
%   Einc = richmond_plane_wave_incident(X, Y, k0, theta_i)
%
%   X, Y    : matrices of coordinates
%   k0      : free-space wavenumber
%   theta_i : incidence angle (from +x axis)
%
%   Einc    : Ez^inc at each (x,y)

kx = k0 * cos(theta_i);
ky = k0 * sin(theta_i);

Einc = exp(-1i * (kx * X + ky * Y));
end
