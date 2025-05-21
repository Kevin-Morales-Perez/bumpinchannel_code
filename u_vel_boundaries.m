function [u_vel] = u_vel_boundaries(u_vel,n_x,n_y)
%Boundary conditions for velocity in X
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index
    
    %outlet boundary condition, zero gradient  (Neumman)
    u_vel(:,n_x-1)=u_vel(:,n_x-2);
    
    %No wall Slip condition (Dirichlet) on bump
    u_vel(n_y-1,le_index:te_index)=-u_vel(n_y-2,105:253);

    %zero shear stress outside the bump (Neumman)
    u_vel(n_y-1,1:le_index)=u_vel(n_y-2,1:le_index);
    u_vel(n_y-1,te_index:n_x-1)=u_vel(n_y-2,te_index:n_x-1);

    %top boundary, zero gradient (Neumman)
    u_vel(1,:)=u_vel(2,:);

end