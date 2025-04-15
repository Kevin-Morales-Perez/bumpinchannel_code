function [d_wall] = near_walld(node_cords)
%Function to calculate nearest distance to the wall
%finding nearest distance wall from cell centroids
x=node_cords(1);
y=node_cords(2);
    %case 1) Node point lies backwars or upwards plate
    if x < 0 || x > 1.5
        if x< 0
            p_wall=[0,0];
            d_wall=sqrt(( x - p_wall(1) )^2 + (y - p_wall(2))^2 );
        else
            p_wall=[1.5,0];
            d_wall=sqrt(( x - p_wall(1) )^2 + (y - p_wall(2))^2 );       
        end
    else
    %case 2) Node is above the plate    
    %use  Newton-Rhapson to find nearest distance to the wall
    %find nearest point using newton raphson
    [x_n,y_n] = Newton_Raphson_nd(x,y);
    %case 2.1) calculate with newton Raphson
    d_walln=sqrt(( x - x_n)^2 + (y - y_n)^2 );
    %case 2.2) just take the y coordinate as the distance 
    d_wall_b= y;
    d_wall=min(d_wall_b,d_walln);%compare 
    end
end