function [tau_xx,tau_xy,tau_yy] = turbulent_stressses_boundaryc(tau_xx,tau_xy,tau_yy,n_x,n_y)
%Boundary conditions for velocity in X
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index

    %######### INLET  ############
    %####### Dirichlet ############
    tau_xx(:,1)=-tau_xx(:,2);
    tau_xy(:,1)=-tau_xy(:,2);
    tau_yy(:,1)=-tau_yy(:,2);

    
    %######### TOP  ############
    %####### Neumman ############
    tau_xx(1,:)=tau_xx(2,:);
    tau_xy(1,:)=tau_xy(2,:);
    tau_yy(1,:)=tau_yy(2,:);

    %######### OUTLET ###########
    %####### Neumman ############
    tau_xx(:,n_x-1)=tau_xx(:,n_x-2);
    tau_xy(:,n_x-1)=tau_xy(:,n_x-2);
    tau_yy(:,n_x-1)=tau_yy(:,n_x-2);

    %####### ON BUMP #############
    %####### Dirichlet ###########
    tau_xx(n_y-1,le_index:te_index)=-tau_xx(n_y-2,le_index:te_index);
    tau_xy(n_y-1,le_index:te_index)=-tau_xy(n_y-2,le_index:te_index);
    tau_yy(n_y-1,le_index:te_index)=-tau_yy(n_y-2,le_index:te_index);

    %#### OUTSIDE THE BUMP ####### 
    %####### Neumman ###########

    %-Upwind
    tau_xx(n_y-1,1:le_index)=-tau_xx(n_y-2,1:le_index);
    tau_xy(n_y-1,1:le_index)=-tau_xy(n_y-2,1:le_index);
    tau_yy(n_y-1,1:le_index)=-tau_yy(n_y-2,1:le_index);
    %Downwind
    tau_xx(n_y-1,te_index:n_x-1)=-tau_xx(n_y-2,te_index:n_x-1);
    tau_xy(n_y-1,te_index:n_x-1)=-tau_xy(n_y-2,te_index:n_x-1);
    tau_yy(n_y-1,te_index:n_x-1)=-tau_yy(n_y-2,te_index:n_x-1);

end