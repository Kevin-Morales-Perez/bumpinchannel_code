%Grid with geometrical inflation layer for bump in channel problem
clear all
close all
clc
%tic

%boundary layer thickness estimation
Re=3e6;%Reynolds Number
x_lenght=1.5;%lenght of the plate
d99=(0.35*x_lenght)/(Re^(1/5)) + 0.05*2;%boundary layer estimation

y_wall =3.8e-4;%Non dimensional wall distance
y_h=2*y_wall;%height of the first cell
inf_layer_N = 15;% Number of layers
%d99=1.0230;%height of the boundary layer

%%% Computig geometrical growth factor
inf_layer_G=0.00001;%Geometrical factor 
y_t=0;%Boundary layer height
while y_t<=d99
    inf_layer_G=inf_layer_G+0.00001;
    y_t=y_h*((1-inf_layer_G^inf_layer_N)/(1-inf_layer_G));
end    

l_x = 50; % length in X
l_y = 5; % length in Y
n_x = 354;  % number of grid points in X
n_y =162;  % number of points in Y

or_plate=[0 0];%origin of the plate X,Y
len_plate=1.5;%len of the plate
p_farfield_l=or_plate(1)- 0.5*(l_x-len_plate);%position of the left farfield
p_farfield_r=or_plate(1)+len_plate+ 0.5*(l_x-len_plate);%position of the right farfield
p_farfield_t=or_plate(2)+l_y;%position of the top farfield

%Creating X steps with refinement at the leading and  trailing edges
%vector with the X grid coordinates
%x=zeros(1,n_x+1);
%y=zeros(1,n_x+1);% for plotting pourposes


%refinement in the leading and trailing edges with inflation layers x
%direction
inf_layer_Nx=30;%simetrical ,for trailling and leading,proposed
x_h=4*y_h;%first  x step  from leading and trailling egdes
inf_layer_G_ed=1.05;%G factor for inflation layer at the edges

%leading edge
x_plate_led=zeros(1,1+2*inf_layer_Nx);%vector for define steps in leading edge
x_plate_led(21)=or_plate(1);%defining grid point for leading edge
sum_layer_x=0;%len of the inflation layer
%loop for define gridpoints in leading edge
for i =0:inf_layer_Nx-1
    sum_layer_x=sum_layer_x-x_h*inf_layer_G_ed^i;
    x_plate_led(inf_layer_Nx-i)=sum_layer_x;
    x_plate_led(inf_layer_Nx+2+i)=-sum_layer_x;
end
    
%trailling edge
x_plate_trl=zeros(1,1+2*inf_layer_Nx);%vector for define steps in trailling edge
x_plate_trl(inf_layer_Nx+1)=or_plate(1)+len_plate;%defining grid point for trailing edge
sum_layer_x=0;% inflation layer position
%loop for define gridpoints in leading edge
for i =0:inf_layer_Nx-1
    sum_layer_x=sum_layer_x-x_h*inf_layer_G_ed^i;
    x_plate_trl(inf_layer_Nx-i)=x_plate_trl(inf_layer_Nx+1)+sum_layer_x;
    x_plate_trl(inf_layer_Nx+2+i)=x_plate_trl(inf_layer_Nx+1)-sum_layer_x;
end

%portion of the grid in x between the 2 inflation layers
%this portion of the plate has the maximum step that have the inflation
%layers at the edges.

dx_int=x_plate_trl(2*inf_layer_Nx+1)-x_plate_trl(2*inf_layer_Nx);%maximum step in inflation layer
len_int_plate=x_plate_trl(1)-x_plate_led(2*inf_layer_Nx+1);%lenght in the intermediate zone of the plate
gp_int_plate= round(len_int_plate/dx_int);%number of grid points in the interior of the plate
dx_int=len_int_plate/gp_int_plate;%recomputing step in the interior of the plate


x_plate_int=zeros(1,gp_int_plate+1);%number of grids in the interior of the plate 
x_plate_int(1)=x_plate_led(2*inf_layer_Nx+1);
for i =2:gp_int_plate+1
    x_plate_int(i)=x_plate_int(i-1)+dx_int;
end

%number of grid points from edge inflation layers to farfields per side
gp_edf=0.5*(n_x - (4*inf_layer_Nx + gp_int_plate));


%distance from farfield to inflation layers
if x_plate_led(1) - p_farfield_l == p_farfield_r- x_plate_trl(2*inf_layer_Nx+1)
    len_ff=x_plate_led(1) - p_farfield_l;
else 
    printf("Error in  the mesh")
end 
 
%size of step from farfields to inflation layers
dx_ff= len_ff/gp_edf;

%grid portion from left farfield to leading inflation layer
x_grid_lft= zeros(1,gp_edf+1);
x_grid_lft(1)=p_farfield_l;
for i=2:gp_edf+1
    x_grid_lft(i)= x_grid_lft(i-1)+dx_ff;
end

%grid portion from trailling inflation layer to right farfield
x_grid_rgt= zeros(1,gp_edf+1);
x_grid_rgt(1)=x_plate_trl(2*inf_layer_Nx+1);
for i=2:gp_edf+1
    x_grid_rgt(i)= x_grid_rgt(i-1)+dx_ff;
end

x_grid_vector=[x_grid_lft x_plate_led(2:2*inf_layer_Nx+1) x_plate_int(2:gp_int_plate+1) x_plate_trl(2:2*inf_layer_Nx) x_grid_rgt];

y=zeros(size(x_grid_vector));
plot(x_grid_vector, y, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

fprintf("Grid points in the plate")
disp(2*inf_layer_Nx+gp_int_plate)

%inflation layer in Y direction

y = zeros(1, n_x+1);

b_1 = int16(0.3/dx);
b_2 = int16(1.2/dx);

for i = b_1:b_2
    y(i+1) = 0.05*(sin(pi*x(i+1)/0.9 - (pi/3.)))^4;%comment this line for a flat mesh
end

d_k = l_y/2;  % Distance from the viscous wall where the mesh will be denser
p_k = 50;     % Percentage concentration for the zone with higher mesh density
n_yp = n_y*(p_k/100); % Number of nodes in the high-density region
n_yf = n_y - n_yp;    % Number of nodes in the free region

dy_p = d_k / n_yp;     % y-step in the high-density part of the mesh
dy_f = (l_y - d_k) / n_yf; % y-step in the remaining part of the mesh

% Matrix that will store the X coordinates of the mesh
X = zeros(n_y+1, n_x+1);
% Storing X coordinates, constant step
for i = 1:n_x+1
    X(:, i) = x(i);
end

% Initial matrix that will store the Y coordinates of the mesh
Y_t = zeros(n_y+1, n_x+1);
for i = 1:n_x+1
    Y_t(1, i) = y(i);
end

% Storing Y coordinates for the high-density nodes
for i = 2:n_yp+1
    for j = 1:n_x+1
        Y_t(i, j) = Y_t(i-1, j) + dy_p;
    end
end

% Storing Y coordinates for the rest of the mesh, dynamic dy_f
dy_f_d = 0;          % dynamic y-step
max_y = max(y);      % maximum value of Y
dy_din = zeros(1, n_x+1); % vector for y-steps in the free part
for i = n_yp+2:n_y+1
    for j = 1:n_x+1
        % dy_f_d = (l_y + max_y - y(j)) / n_yf;
        dy_f_d = (l_y - (y(j) + d_k)) / n_yf;
        Y_t(i, j) = Y_t(i-1, j) + dy_f_d;
        dy_din(j) = dy_f_d / dx;
    end
end

Y = zeros(size(Y_t));

k = 0;
for i = 1:n_y+1
    for j = 1:n_x+1
        Y(i, j) = Y_t(n_y+1-k, j);
    end
    k = k + 1;
end

clear Y_t;

z = zeros(n_y+1, n_x+1);

figure(1)
plot(x, y)

% Flatten the matrices into vectors
x_coords = X(:);
y_coords = Y(:);
%{
% Create a scatter plot
figure(2)
scatter(x_coords, y_coords, 10, 'filled');  % The '10' specifies the marker size
% Label the axes
xlabel('X Coordinates');
ylabel('Y Coordinates');
% Title for the plot
title('Scatter Plot of Unstructured Mesh Points');
% Set axis equal for better representation of the points
axis equal;
% Optionally, add grid and adjust other plot settings
grid on;
%}
figure(3)
mesh(X, Y, z)
title('Mesh')
axis equal
