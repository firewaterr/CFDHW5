function Psi = iterateStreamfunctionField(Psi_init, relaxationFactor, Nx, Ny, Lx, Ly, Omega, tolerance)
    % Initialize variables
    NormError = inf;
    iteration = 0;
    Psi = Psi_init;
    delta = Lx / Nx; % Grid spacing
    % Main loop
    while NormError > tolerance
        iteration = iteration + 1; % Update iteration counter
        Psi_old = Psi;
        Psi_relax = Psi;
        % Forward update
        for i = 2:Nx-1
            for j = 2:Ny-1
                Psi_relax(i,j) = (Psi_relax(i+1,j) + Psi_relax(i-1,j) + Psi_relax(i,j+1) + Psi_relax(i,j-1) + Omega(i,j) * delta^2) / 4;
            end
        end

        % Backward update
        for i = Nx-1:-1:2
            for j = Ny-1:-1:2
                Psi_relax(i,j) = (Psi_relax(i+1,j) + Psi_relax(i-1,j) + Psi_relax(i,j+1) + Psi_relax(i,j-1) + Omega(i,j) * delta^2) / 4;
            end
        end

        % Relaxation method
        for i = 2:Nx-1
            for j = 2:Ny-1
                Psi(i,j) = (1 - relaxationFactor) * Psi(i,j) + relaxationFactor * Psi_relax(i,j);
            end
        end

        % Calculate error
        error = (Psi - Psi_old).^2;
        NormError = sqrt(mean(error(:))); % Calculate the error norm
    end
end