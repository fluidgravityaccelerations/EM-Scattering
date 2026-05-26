function [Xc, Yc, dx, dy] = richmond_make_grid(Lx, Ly, Nx, Ny)
%RICHMOND_MAKE_GRID Uniform Cartesian grid of cell centers.
%
%   [Xc, Yc, dx, dy] = richmond_make_grid(Lx, Ly, Nx, Ny)
%
%   Lx, Ly : domain size in x and y
%   Nx, Ny : number of cells in x and y
%
%   Xc, Yc : matrices of cell-center coordinates
%   dx, dy : cell sizes

dx = Lx / Nx;
dy = Ly / Ny;

xv = linspace(-Lx/2 + dx/2, Lx/2 - dx/2, Nx);
yv = linspace(-Ly/2 + dy/2, Ly/2 - dy/2, Ny);

[Xc, Yc] = meshgrid(xv, yv);
end
