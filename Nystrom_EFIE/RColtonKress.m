function out = RColtonKress(N, t, tj)
%RColtonKress Colton–Kress trigonometric weight for the periodic log kernel.
%
%   out = RColtonKress(N, t, tj)
%
%   N   : number of intervals (2N nodes on [0, 2*pi))
%   t   : target angle
%   tj  : source angle

out = 0;
for m = 1:N-1
    out = out + (1/m)*cos(m*(t - tj));
end
out = -(2*pi/N)*out - pi/(N^2)*cos(N*(t - tj));
end
