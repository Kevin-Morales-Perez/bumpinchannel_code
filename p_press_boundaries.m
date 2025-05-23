function [p_press] = p_press_boundaries(p_press,n_x,n_y,y_ghost_cells,...
    x_wall_cells,y_wall_cells)
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
    p_press(:,n_x-1)=1.02828*p_press(:,n_x-2);% To be reviewed 

    %#### OUTSIDE THE BUMP ####### 
    %####### Neumman ###########
    %inlet and outlow flow 
    p_press(n_y-1,:)=p_press(n_y-2,:);
 

    %############   ON BUMP  ##################
    %## Neumman(Normal to the wall) ###########
    %NEW CODE TO ENHACE NEUMMAN BOUNDARY CONDITION  ON BUMP
    %GET NODE COORDINATES OF CELLS THAT ARE ON THE WALL
    %GET NODE COORDINATES OF BOTTOM GHOST CELLS (BENEATH THE WALL)
    for i=le_index:te_index
        %determine sign of  SLOPE in ghost cells
        slope=y_ghost_cells(i+1)-y_ghost_cells(i);
        if slope >0
            theta_p=atan(slope);
            delta_y=y_wall_cells(i+1)-y_ghost_cells(i+1);
            dist_a=delta_y*sin(theta_p);
            distance_p=sqrt((x_wall_cells(i+1)-x_wall_cells(i))^2+...
                (y_wall_cells(i+1)-y_wall_cells(i))^2);
            diff_press=p_press(n_y-2,i+1)-p_press(n_y-2,i);
            p_inter=p_press(n_y-2,i) + (diff_press/distance_p)*(distance_p-...
                dist_a);
            p_press(n_y-1,i+1)=p_inter;
       
        elseif slope<0
            theta_p=atan(slope);
            delta_y=y_wall_cells(i)-y_ghost_cells(i);
            dist_a=delta_y*sin(theta_p);
            distance_p=sqrt((x_wall_cells(i+1)-x_wall_cells(i))^2+...
                (y_wall_cells(i+1)-y_wall_cells(i))^2);
            diff_press=p_press(n_y-2,i+1)-p_press(n_y-2,i);
            p_inter=p_press(n_y-2,i) + (diff_press/distance_p)*dist_a;
            p_press(n_y-1,i+1)=p_inter;
        else
            p_press(n_y-1,i+1)=p_press(n_y-2,i+1);
        end
    end
end

