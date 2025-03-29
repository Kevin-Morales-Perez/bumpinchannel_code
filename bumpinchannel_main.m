
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bump in channel problem
%Written by Kevin Morales 
%Instituto Politécnico Nacional , Aeronautical Engineering
clear all
close all 
clc
tic

%________________________Domain____________________________
%Load mesh
%X for x an Y for Y coordinates of the points of the mesh
load("mesh_bumpchannel2.mat");
%load("mesh_bumpchannel353_161.mat");

[n_y,n_x] =size(X);%number of cells in  X and Y
dx=X(1,2)-X(1,1);% X step 
l_y = max(Y(:,1));%Lenght in Y

%_________________geometrical properties___________________
x_cent= zeros(n_y-1,n_x-1); % matrix which contains x centroids of cells  
y_cent= zeros(n_y-1,n_x-1); %matrix which contains y centroids of cells 
len_faces=zeros(n_y-1,n_x-1,4); %matrix which contains lenght of the 4 faces of each cell
trig_cells=zeros(n_y-1,n_x-1,2,3);%matrix which contains trig funcions for angles in faces 2) and 4)
face_centers=zeros(n_y-1,n_x-1,4,2);%matrix which contains the coordinates of points at the centres of each face
cell_volumes=zeros(n_y-1,n_x-1);%matrix which contains the volumen of each cell


%x3=[0 0 4 4]
%y3=[0 4 5 2]
%cell_collocated_node(x3,y3)

for i=1:n_y-1
    for j=1:n_x-1
        x1=[X(i+1,j) X(i,j) X(i,j+1) X(i+1,j+1)];
        y1=[Y(i+1,j) Y(i,j) Y(i,j+1) Y(i+1,j+1)];
        [x_cent(i,j),y_cent(i,j),len_faces(i,j,:),trig_cells(i,j,:,:),face_centers(i,j,:,:),cell_volumes(i,j)]=cell_collocated_node(x1,y1);
    end
end

%________________________Constants________________________
vel_ini=9e-5;  %29e-4; %Velocidad en la entrada 
rho =1.2; %Density (Kg/m3)
nu =0.0000174; %Viscosity (Kg/(m*s))
reynolds = vel_ini*rho*l_y/nu;%Reynolds number 
flux_mass=l_y*vel_ini;%mass flux
disp(reynolds)

%_______________Velocity and presure fields_______________
u_vel=zeros(n_y+1,n_x+1); %Velocity in X axis
v_vel=zeros(n_y+1,n_x+1); %Velocity in Y axis
press=ones(n_y+1,n_x+1);  %Pressure

%_______________underrelaxation_factors____________________
alpha_p=0.5;
alpha_v=0.5;

%________________________err________________________________
b_p =zeros(n_y+1,n_x+1);%mass balance in each cell
err=1;%error as residual of continuity
err_req=1e6; %error required
err_vec=[];%vector to storage the errors
iter=0;%iterations



%{
%_____________ main iteration  ______________
while err > err_req
    if err > 1e12
        fprintf("Unstable solution")
        break
    end
    %Ecuacion de momento en X, Y , y corrección de presiones
    for i=1:n_y
        for j =1:n_x-1
            %xmomentum
        end
    end

   %x boundary conditions
   for i=1:n_y
        for j =1:n_x-1
            %ymomentum
        end
    end
    %y_boundary
    %pressure correction
    %calculating err
end
%}

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
