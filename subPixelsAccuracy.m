function [p_y,p_x] = subPixelsAccuracy( matrix, p_y_in, p_x_in, win )

[height_dev,width_dev] = size(matrix);

points = zeros(win*win,1);
coeff_matrix = zeros(win*win,6);
nb_points = size( p_x_in,2 );
p_y = p_y_in;
p_x = p_x_in;

for i = 1:nb_points
    %%% Compute the neighbourhood
    neigh_ind = 1;
    for j = p_y_in(i) - (win-1)/2:p_y_in(i) + (win-1)/2
        for k = p_x_in(i) - (win-1)/2:p_x_in(i) + (win-1)/2
            if (j > 1)&&(j < height_dev)&&(k > 1)&&(k < width_dev)
                points(neigh_ind) = matrix(j,k);
                neigh_ind = neigh_ind + 1;
            end
        end
    end
    coeff_matrix(1,:) = [1,1,-1,-1,1,1];
    coeff_matrix(2,:) = [0,1,0,-1,0,1];
    coeff_matrix(3,:) = [1,1,1,-1,-1,1];
    
    coeff_matrix(4,:) =[1,0,-1,0,0,1];
    coeff_matrix(5,:) =[0,0,0,0,0,1];
    coeff_matrix(6,:)=[1,0,1,0,0,1];
    
    coeff_matrix(7,:)=[1,1,-1,1,-1,1];
    coeff_matrix(8,:)=[0,1,0,1,0,1];
    coeff_matrix(9,:)=[1,1,1,1,1,1];
    
    
    %%% Compute the coefficient
    coeff = ((coeff_matrix' * coeff_matrix)^-1) * coeff_matrix' * points;
    
    %%% Compute the offset
    dev_1_mat = [2*coeff(1), coeff(5) ; coeff(5), 2*coeff(2)];
    dev_2_mat = [-coeff(3) ; -coeff(4)];
    offset = (dev_1_mat'*dev_1_mat)^-1*dev_1_mat*dev_2_mat;
    p_x(i) = p_x_in(i) + offset(1);
    p_y(i) = p_y_in(i) + offset(2);
end
