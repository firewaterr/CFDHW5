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
Nt = 200;
dt = T / Nt; % Time step
mu = 0.001; % Viscosity
relaxationFactor = 1.95; % Relaxation factor for streamfunction iteration
%% Initial Conditions
Psi     = zeros(Nx+2, Ny+2); % Streamfunction
Omega   = zeros(Nx+2, Ny+2); % Vorticity
U       = zeros(Nx+2, Ny+2); % x-velocity
V       = zeros(Nx+2, Ny+2); % y-velocity
t       = linspace(0, T, Nt); % Time vector
x       = linspace(Lx, 0, Nx); % x-coordinates
y       = linspace(0, Ly, Ny); % y-coordinates

%% Boundary Conditions
vTop    = sin(pi .* x).^2;  % Top wall velocity
vBottom = zeros(1, Nx);     % Bottom wall velocity
vLeft   = zeros(Ny, 1);     % Left wall velocity
vRight  = zeros(Ny, 1);     % Right wall velocity

%% Main Loop
Loopout = 1e-4;
iteration = 0;
maxIterations = 100000;
deltaOmega = inf;
while (deltaOmega > Loopout) && (iteration < maxIterations)
    iteration = iteration + 1;

    %% Step 1: Calculate Vorticity for Inner Points
    Omega_old = Omega;
    for i = 2:Nx+1
        for j = 2:Ny+1
            Omega(i,j) = Omega_old(i,j) + dt * (...
                (Psi(i,j+1) - Psi(i,j-1)) / (2*dy) * (Omega_old(i+1,j) - Omega_old(i-1,j)) / (2*dx) - ...
                (Psi(i+1,j) - Psi(i-1,j)) / (2*dx) * (Omega_old(i,j+1) - Omega_old(i,j-1)) / (2*dy) + ...
                mu * ((Omega_old(i+1,j) - 2*Omega_old(i,j) + Omega_old(i-1,j)) / dx^2 + ...
                   (Omega_old(i,j+1) - 2*Omega_old(i,j) + Omega_old(i,j-1)) / dy^2));
        end
    end
    %% Step 2: Calculate Streamfunction
    Psi = iterateStreamfunctionField(Psi, relaxationFactor, Nx+2, Ny+2, Lx, Ly, Omega, 1e-6);

    %% Step 3: Calculate Vorticity for Boubndary Points
    for i = 2:Ny+1
        Omega(1, i)   = -2 * (Psi(2, i) - Psi(1, i) + vLeft(i-1) * dx) / dx^2; % Left wall
        Omega(end, i) = -2 * (Psi(end-1, i) - Psi(end, i) + vRight(i-1) * dx) / dx^2; % Right wall
    end
    for j = 2:Nx+1
        Omega(j, 1) = -2 * (Psi(j, 2) - Psi(j, 1) + vBottom(j-1) * dy) / dy^2; % Bottom wall
        Omega(j, end)   = -2 * (Psi(j, end-1) - Psi(j, end) + vTop(j-1) * dy) / dy^2; % Top wall
    end
    error = (Omega - Omega_old).^2;
    deltaOmega = sqrt(mean(error(:)));
end

%% Step 4: Calculate Velocity
for i = 2:Nx+1
    for j = 2:Ny+1
        U(i,j) = (Psi(i,j+1) - Psi(i,j-1)) / (2*dy);
        V(i,j) = (Psi(i+1,j) - Psi(i-1,j)) / (2*dx);
    end
end

%% Center of the Vorticity
[minValue, minIndex] = min(Psi(:));
[rowIndex, colIndex] = ind2sub(size(Psi), minIndex);
Centerx = (Nx+2-rowIndex)/(Nx+2);
Centery = colIndex/(Ny+2);
disp(['最小值: ', num2str(minValue)]);
disp(['横坐标 (行): ', num2str(Centerx)]);
disp(['纵坐标 (列): ', num2str(Centery)]);
[maxValue, maxIndex] = max(Psi(:));
[rowIndex, colIndex] = ind2sub(size(Psi), maxIndex);
Centerx = (Nx+2-rowIndex)/(Nx+2);
Centery = colIndex/(Ny+2);
disp(['最大值: ', num2str(maxValue)]);
disp(['横坐标 (行): ', num2str(Centerx)]);
disp(['纵坐标 (列): ', num2str(Centery)]);

%% Plot the Results
figure;
contourf(x ,y, Psi(2:Nx+1,2:Ny+1)',20,'EdgeColor','none');
%view(2);
colorbar;
colormap("parula");
xlabel('x');
ylabel('y');
hold on;
l = streamslice(x, y, U(2:Ny+1,2:Nx+1)', V(2:Ny+1,2:Nx+1)',5);
set(l,'LineWidth',1)
set(l,'Color','black');
title('Streamlines of Velocity Field');
axis equal tight;
U_middlex = U(round(Nx/2), 2:Ny+1);
V_middlex = V(round(Nx/2), 2:Ny+1);
U_middley = U(2:Nx+1, round(Ny/2));
V_middley = V(2:Nx+1, round(Ny/2));
figure;
plot(y, U_middlex, 'r', 'LineWidth', 1.5);
hold on;
plot(y, V_middlex, 'b', 'LineWidth', 1.5);
xlabel('y');
ylabel('Velocity');
title('Velocity Profile along the Middle x of the Cavity');
legend('U', 'V');
figure;
plot(x, U_middley, 'r', 'LineWidth', 1.5);
hold on;
plot(x, V_middley, 'b', 'LineWidth', 1.5);
xlabel('x');
ylabel('Velocity');
title('Velocity Profile along the Middle y of the Cavity');
legend('U', 'V');