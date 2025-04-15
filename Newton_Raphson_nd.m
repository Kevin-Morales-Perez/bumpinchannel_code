function [x_sol,f_x] = Newton_Raphson_nd(px,py)
    %Newton - Raphson method to obtain calculate nearest distance to the wall
    %Newton -Raphson aproximation to solve dot product equals to zero, between
    %the distance vector from node to nearest point of the wall , and unitary
    %tangent vector in this point , which is the derivative at that point
    %valid range from 0.3 to 1.2 
    x_sol=0.3;%x initial
    g_x=1;
    dA_dx=pi*0.9;
    k=0.2*pi/0.9;
    if py>2
        err=0.8;
    elseif py>1
        err=0.2;
    else
        err=1e-10;
    end
    
    while (g_x > err)
        A=((pi*x_sol)/0.9) - pi/3;%argument of main function that describes the wall
        f_x =0.05*(sin(A))^4;% main function that describes the wall
        %---------------------Below just derivatives -------------------------
        df_dx=k*((sin(A))^3)*cos(A);
        df2_dx2=k*dA_dx*(3*(((sin(A))^2)*((cos(A))^2)) - (sin(A))^4);
        g_x=(px-x_sol)+(py-f_x)*df_dx;
        dg_dx=-1 +(py-f_x)*df2_dx2 - (df_dx^2);
        x_sol=x_sol - g_x/dg_dx;%main Newton Raphson
    end
end