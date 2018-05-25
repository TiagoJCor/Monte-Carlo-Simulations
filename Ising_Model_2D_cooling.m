% Ising_model_2D_cooling.m

% This script performs several simulations of a Monte Carlo 2D Ising Model 
% scheme at different temperatures in order to find the critical
% temperature at which a state change occurs, and the value of the global
% system magnetization.

% Work by Tiago Correia 
% Winter 2015

% Size of the grid (note that each grid has grid_size^2 elements.
grid_side = 4;      
% Number of steps each simulation is allowed to run for.
nsteps = 6000;      

% Define distribution.
f = @ (E,Ts) exp(-E/Ts);
% Number of simulations to run.
n = 200;             

for p = 1:n
    
    energy = zeros(nsteps,1);
    T(p) = rand()*5 + 1e-5;
    grid = ones(grid_side, grid_side);
    
    % Initialize lattice
    for i = 1:grid_side  
        for j = 1:grid_side
            % Populate lattice with +1/-1 (random).
            grid(i,j) = grid(i,j) - 2 * fix(2*rand);
        end
    end
    
    % Calculate initial energy.
    energy(1) = calc_energ2d(grid);
    accept = 0;
    for t = 2:nsteps
        % Choose row and column of index to flip (randomly).
        flip_i = randi(grid_side);
        flip_j = randi(grid_side);
        
        % Propose new state.
        grid_new = grid;
        grid_new(flip_i, flip_j) = -grid_new(flip_i, flip_j);
        % Calculate new energy.
        H_new = single_calc_energ2d(grid_new, flip_i, flip_j);
        % Record old energy.
        H_old = single_calc_energ2d(grid, flip_i, flip_j);
        
        h = min(1,f((H_new - H_old), T(p)));
        
        % Accept state change, update grid and energy.
        if (rand < h)      
            grid = grid_new;
            energy(t) = energy(t-1) + (H_new - H_old);
            accept = accept+1;
        else
            % Don't accept, grid and energy stay the same.
            energy(t) = energy(t-1);  
        end

    end

    % Calculate average energy and magnetization for this given
    % temperature.
    avg_E(p) = energy(nsteps)/grid_side^2;
    avg_M(p) = sum(sum(grid))/grid_side^2;
end

%=========================================================================%

% Plots

figure(2)
plot(T, avg_M, 'or')
ylabel('Average Magnetization')
xlabel('Temperature')
title([num2str(grid_side) ' by ' num2str(grid_side) ' lattice, ',...
    num2str(grid_side) ' data points, ' num2str(nsteps) ' steps per simulation'])

figure(3)
plot(T, avg_E, 'ob')
ylabel('Average Energy')
xlabel('Temperature')
title([num2str(grid_side) ' by ' num2str(grid_side) ' lattice, ',...
    num2str(grid_side) ' data points, ' num2str(nsteps) ' steps per simulation'])