% CFDHW5                                        %
% Wang Yufeng                                   %
% 2200011013                                    %
% 2025.05.09                                    %
% Lid-driven cavity flow                        %
% 2D incompressible Navier-Stokes equations     %   
% Vorticity-Streamfunction Method               %
% Boundary Conditions                           %
% Top wall:     U = (sin(\pi * x))^2,   V = 0   %
% Bottom wall:  U = 0,                  V = 0   %
% Left wall:    U = 0,                  V = 0   %
% Right wall:   U = 0,                  V = 0   %

%% Code Begin
clear; 
clc; 
close all;

%% Parameters
Nx = 100;
Ny = 100;
Lx = 1;
Ly = 1;
dx = Lx / Nx;
dy = Ly / Ny;
T = 1; % Total time
Nt = 100;
dt = T / Nt; % Time step
mu = 0.001; % Viscosity
relaxationFactor = 1.9; % Relaxation factor for streamfunction iteration
%% Initial Conditions
Psi     = zeros(Ny+2, Nx+2); % Streamfunction
Omega   = zeros(Ny+2, Nx+2); % Vorticity
U       = zeros(Ny+2, Nx+2); % x-velocity
V       = zeros(Ny+2, Nx+2); % y-velocity
t       = linspace(0, T, Nt); % Time vector
x       = linspace(0, Lx, Nx+2); % x-coordinates
y       = linspace(0, Ly, Ny+2); % y-coordinates

%% Boundary Conditions
vTop    = sin(pi .* x).^2;  % Top wall velocity
vBottom = zeros(1, Nx);     % Bottom wall velocity
vLeft   = zeros(Ny, 1);     % Left wall velocity
vRight  = zeros(Ny, 1);     % Right wall velocity

%% Main Loop
for k = 1:Nt
    %% Step 1: Calculate Vorticity for Inner Points
    Omega_old = Omega;
    for i = 2:Ny+1
        for j = 2:Nx+1
            Omega(i,j) = Omega_old(i,j) + dt * (...
                (Psi(i+1,j) - Psi(i-1,j)) / (2*dx) * (Omega_old(i,j+1) - Omega_old(i,j-1)) / (2*dy) - ...
                (Psi(i,j+1) - Psi(i,j-1)) / (2*dy) * (Omega_old(i+1,j) - Omega_old(i-1,j)) / (2*dx) + ...
                mu * ((Omega_old(i+1,j) - 2*Omega_old(i,j) + Omega_old(i-1,j)) / dx^2 + ...
                   (Omega_old(i,j+1) - 2*Omega_old(i,j) + Omega_old(i,j-1)) / dy^2));
        end
    end

    %% Step 2: Calculate Streamfunction
    Psi = iterateStreamfunctionField(Psi, relaxationFactor, Nx+2, Ny+2, Lx, Ly, Omega, 1e-6);

    %% Step 3: Calculate Vorticity for Boubndary Points

    for i = 2:Ny+1
        Omega(i, 1)   = -2 * (Psi(i, 2) - Psi(i, 1) + vLeft(i-1) * dx) / (dx)^2; % Left wall
        Omega(i, end) = -2 * (Psi(i, end-1) - Psi(i, end) + vRight(i-1) * dx) / (dx)^2; % Right wall
    end
    for j = 2:Nx+1
        Omega(1, j)   = -2 * (Psi(2, j) - Psi(1, j) + vBottom(j-1) * dy) / (dy)^2; % Bottom wall
        Omega(end, j) = -2 * (Psi(end-1, j) - Psi(end, j) + vTop(j-1) * dy) / (dy)^2; % Top wall
    end

end
%% Plot the Results

figure;
contourf(x, y, Omega, 20);
view(2);
colorbar;
xlabel('x');
ylabel('y');
title('Vorticity Distribution');
