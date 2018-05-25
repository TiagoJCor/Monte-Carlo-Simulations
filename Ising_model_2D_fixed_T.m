% Ising_model_2D_fixed_T.m

% This script performs a single simulation of a Monte Carlo 2D Ising Model 
% scheme at a given temperature and plots the resulting system.

% Work by Tiago Correia 
% Winter 2015

% Number of steps the simulation is allowed to run for.
nsteps = 10^2;      
% Fix temperature.
beta = 100;         
% Initialize magnetization.
M = 0;              

d = zeros(nsteps,1);
H = zeros(nsteps,1);

% Size of the grid (note that each grid has grid_size^2 elements.
grid_side = 4;      
grid = ones(grid_side,grid_side); 

for i = 1:grid_side  
    for j = 1:grid_side
        % Populate lattice with +1/-1 (random).
        grid(i,j) = grid(i,j) - 2*fix(2*rand);
    end
end

% Plots initial configuration.

figure(1)
title('Initial Configuration')
hold on
axis off
for i = 1:grid_side
    for j = 1:grid_side
        if(grid(i,j) == 1)  
            % Up mag. moment.
            plot(i,j,'+r')
        else
            % Down mag. moment.
            plot(i,j,'.k')  
        end
    end
end

% Calculate initial energy.
H(1) = calc_energ2d(grid);
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
    
    h = min(1,exp(-beta*(H_new - H_old)));
    
    % Accept state change, update grid and energy.
    if (rand < h)      
        grid = grid_new;
        H(t) = H(t-1) + (H_new - H_old);
        accept = accept+1;
    else
        % Don't accept, grid and energy stay the same.
        H(t) = H(t-1);  
    end
    Mx = sum(sum(grid));
    % To obtain avg mag. mom.
    M = M + Mx;        
    d(t) = Mx;
end

figure(2)
histogram(d)
title('Histogram of Magnetic Moment Distribution')
ylabel('Frequency')
xlabel('Magnetic Moment')
figure(3)
hold on
axis off

% Plots final configuration of the system.
for i = 1:grid_side
    for j = 1:grid_side
        if(grid(i,j) == 1)  
            % Up mag. moment.
            plot(i,j,'+r')
        else
            % Down mag. moment.
            plot(i,j,'.k')  
        end
    end
end
title('Final Configuration')

