function [u_vel] = u_vel_boundaries(u_vel,n_x,n_y,vel_inlet)
%Boundary conditions for velocity in X
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index

    %######### INLET  ############
    %####### Dirichlet ############
    u_vel(:,1)=2*vel_inlet - u_vel(:,2);


    %######### TOP  ############
    %####### Neumman ############
    u_vel(1,:)=u_vel(2,:);

    %######### OUTLET ###########
    %####### Neumman ############
    u_vel(:,n_x-1)=u_vel(:,n_x-2);

    %####### ON BUMP #############
    %####### Dirichlet ###########
    u_vel(n_y-1,le_index:te_index)=-u_vel(n_y-2,105:253);


    %#### OUTSIDE THE BUMP ####### 
    %####### Neumman ###########

    %-Upwind
    u_vel(n_y-1,1:le_index)=u_vel(n_y-2,1:le_index);
    %-Downwind
    u_vel(n_y-1,te_index:n_x-1)=u_vel(n_y-2,te_index:n_x-1);


end