function [x, y, xp, yp, tp, Deltal] = discretize_cylinder(a, NN)
%DISCRETIZE_CYLINDER  Discretization of a circular PEC cylinder for 2D MoM.
%
%   [x, y, xp, yp, tp, Deltal] = discretize_cylinder(a, NN)
%
%   Inputs:
%       a   - Cylinder radius
%       NN  - Number of discretization segments (panels)
%
%   Outputs:
%       x, y    - Coordinates of the segment endpoints (size NN)
%       xp, yp  - Coordinates of the matching (collocation) points (size NN)
%       tp      - Angular position of the matching points (size NN)
%       Deltal  - Effective segment length used in MoM integration
%
%   Notes:
%   - The cylinder is discretized into NN straight segments.
%   - Matching points are located at the midpoints of each segment.
%   - The segment length Δl is computed using the accurate geometric
%     expression for midpoint matching:
%
%         Δl = 2 * a * tan(π / (NN + 1))
%
%     This accounts for the fact that the matching points lie halfway
%     between the angular nodes, providing a more precise approximation
%     of the effective integration length.

    % Angular nodes (NN+1 points, last excluded)
    t = linspace(0, 2*pi - 1e-9, NN + 1);

    % Segment endpoints
    x = a * cos(t(1:NN));
    y = a * sin(t(1:NN));

    % Matching points (midpoints of each segment)
    xp(1:NN-1) = (x(2:NN) + x(1:NN-1)) / 2;
    yp(1:NN-1) = (y(2:NN) + y(1:NN-1)) / 2;
    tp(1:NN-1) = (t(2:NN) + t(1:NN-1)) / 2;

    % Last matching point (wrap-around)
    xp(NN) = (x(NN) + x(1)) / 2;
    yp(NN) = (y(NN) + y(1)) / 2;
    tp(NN) = mod((t(1) + (t(NN) - 2*pi)) / 2, 2*pi);

    % Effective segment length
    Deltal = 2 * a * tan(pi / (NN + 1));

end
