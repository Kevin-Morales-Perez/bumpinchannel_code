function [p_press] = p_press_boundaries(p_press,n_x,n_y)
%Boundary conditions for pressure
%Bump in channel problem 

    le_index=105;%Leading edge index 
    te_index=253;%Trailing edge index
   
    %######### INLET  ############
    %#######  INFLOW ############
    p_press(:,1)=(1/1.02828)*p_press(:,2);

    %######### TOP  ############
    %####### Neumman ############
    p_press(1,:)=p_press(2,:);


    %######### OUTLET ###########
    %#######  OUTFLOW  ############
    %p_press(:,n_x-1)=p_press(:,n_x-2);% To be reviewed 

    %############   ON BUMP  ##################
    %## Neumman(Normal to the wall) ###########
    p_press(n_y-1, le_index:te_index) = p_press(n_y-2, le_index:te_index);%to be reviewed but try

    %#### OUTSIDE THE BUMP ####### 
    %####### Neumman ###########
    %inlet flow %to be reviewed
    p_press(n_y-1,1:le_index)=p_press(n_y-2,1:le_index);
    %outlow flow%to be reviewed
    p_press(n_y-1,te_index:n_x-1)=p_press(n_y-2,te_index:n_x-1);

end

