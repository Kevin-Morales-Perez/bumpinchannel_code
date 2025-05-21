function [p_press] = p_press_boundaries(p_press,n_x,n_y)
%Boundary conditions for pressure
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index
   
    %inlet flow %to be reviewed 
    %outlow flow%to be reviewed

    %bump(neuman)
    p_press(n_y-1, node_le:node_te) = p_press(n_y-2, node_le:node_te);%to be reviewed but try
    %outside the bump (Neumman)
    p_press(n_y-1,1:le_index)=p_press(n_y-2,1:le_index);
    p_press(n_y-1,te_index:n_x-1)=p_press(n_y-2,te_index:n_x-1);

    %Top (Neumman)
    p_press(1,:)=p_press(2,:);


end

