function [nu_tilde] = nu_tilde_boundary_conditions(nu_tilde,n_x,n_y)
%Boundary conditions for velocity in Y
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index

    %######### INLET  ############
    %####### Dirichlet ############
    nu_tilde(:,1)=- nu_tilde(:,2);


    %######### TOP  ############
    %####### Neumman ############
    nu_tilde(1,:)=nu_tilde(2,:);

    %######### OUTLET ###########
    %####### Neumman ############
    nu_tilde(:,n_x-1)=nu_tilde(:,n_x-2);

    %####### ON BUMP #############
    %####### Dirichlet ###########

    nu_tilde(n_y-1,le_index:te_index)=-nu_tilde(n_y-2,105:253);


    %#### OUTSIDE THE BUMP ####### 
    %####### Neumman ###########

    %-Upwind
    nu_tilde(n_y-1,1:le_index)=nu_tilde(n_y-2,1:le_index);
    %-Downwind
    nu_tilde(n_y-1,te_index:n_x-1)=nu_tilde(n_y-2,te_index:n_x-1);

end