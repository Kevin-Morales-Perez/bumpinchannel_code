function [v_vel] = v_vel_boundaries(v_vel,n_x,n_y)
%Boundary conditions for velocity in Y
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index
    
    %outlet boundary condition, zero gradient  (Neumman)
    v_vel(:,n_x-1)=v_vel(:,n_x-2);
    
    %No mass flux at the wall (Dirichlet)
    v_vel(n_y-1,le_index:te_index)=-v_vel(n_y-2,105:253);

    %zero gradient outside the bump(Neumman)
    v_vel(n_y-1,1:le_index)=v_vel(n_y-2,1:le_index);
    v_vel(n_y-1,te_index:n_x-1)=v_vel(n_y-2,te_index:n_x-1);

    %top boundary, zero gradient (Neumman)
    v_vel(1,:)=v_vel(2,:);

end