function H = calc_energ2d(grid)
% Calculates magnetization of a given lattice 'grid'.

% Tiago Correia
% Winter 2015

grid_side = length(grid);
H = 0;

for i = 1:grid_side  
    for j = 1:grid_side
        % With periodic boundary conditions.
        if(i == 1)
            a = grid_side;
        else
            a = i-1; 
        end   
        if(i == grid_side)
            b = 1;
        else
            b = i+1;
        end
        
        if(j == 1)
            c = grid_side;
        else
            c = j-1; 
        end
        
        if(j == grid_side)
            d = 1;
        else
            d = j+1;
        end
        H = H - 0.5*(grid(a,j) + grid(b,j) + grid(i,c) + grid(i,d))*grid(i,j); 
    end
end