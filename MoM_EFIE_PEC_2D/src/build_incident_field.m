function Ei = build_incident_field(xp, beta, E0)
%BUILD_INCIDENT_FIELD  Computes the incident plane wave at matching points.
%
%   Ei = build_incident_field(xp, beta, E0)
%
%   Inputs:
%       xp   - x-coordinates of the matching points (size NN)
%       beta - Wavenumber
%       E0   - Amplitude of the incident plane wave
%
%   Output:
%       Ei   - Incident electric field evaluated at the matching points
%
%   Notes:
%   - The incident field is assumed to be a plane wave propagating along +x:
%
%         E_inc(x) = E0 * exp(-j * beta * x)
%
%   - This function can be easily extended to support:
%         * arbitrary incidence angle
%         * cylindrical waves
%         * TE/TM polarization variants
%

    Ei = E0 * exp(-1i * beta * xp(:));

end
