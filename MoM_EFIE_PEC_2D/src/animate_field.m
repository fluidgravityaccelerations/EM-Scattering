function animate_field(Xe, Ye, E, f, ae, be)
%ANIMATE_FIELD  Animates the real part of the total field over time.
%
%   animate_field(Xe, Ye, E, f, ae, be)
%
%   Inputs:
%       Xe, Ye  - Evaluation grid coordinates (size Nx x Ny)
%       E       - Complex total field on the grid (Nx x Ny)
%       f       - Operating frequency
%       ae, be  - Plot extents along x and y axes
%
%   Notes:
%   - The animation shows Re{ E * exp(j ω t) } over one period.
%   - The animation stops when the user presses any key.
%   - The time sampling is chosen to produce a smooth animation.
%

    % Angular frequency
    omega = 2 * pi * f;

    % Time sampling for one period
    dt = 1 / (8 * f);
    omegat = linspace(0, 2*pi, ceil((1 / f) / dt));
    k = 1;

    % Create figure and set keypress callback
    figure;
    set(gcf, 'KeyPressFcn', @(src,event) setappdata(src,'stop',true));
    setappdata(gcf,'stop',false);

    % Animation loop
    while ~getappdata(gcf,'stop')
        surf(Xe, Ye, 2 * real(E * exp(1i * omegat(k))));
        view(2), axis([-ae ae -be be]), axis equal, shading interp, colorbar
        title('Press any key to stop animation');
        drawnow;

        % Advance frame index
        k = mod(k, length(omegat)) + 1;
    end

end
