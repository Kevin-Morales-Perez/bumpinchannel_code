%Bump in channel problem
%Written by Kevin Morales 
%Instituto Politécnico Nacional , Aeronautical Engineering
clear all
%close all 
clc
tic
%________________________Domain____________________________
%Load mesh
%X for x an Y for Y coordinates of the points of the mesh
load("mesh_bumpchannel.mat");
%load("mesh_bumpchannel0_10test.mat");
%load("mesh_bumpchannel353_161.mat");

[n_y,n_x] =size(X);%size of X and Y matrices
l_y = max(Y(:,1));%Lenght in Y
%_________________geometrical properties___________________
cell_centroids=zeros(n_y-1,n_x-1,2);%matrix with the centroids [x,y];
%reshape(cell_centroids(5,5,:),[1,2])
len_faces=zeros(n_y-1,n_x-1,4); %matrix which contains lenght of the 4 faces of each cell, 
%reshape(len_faces(5,5,:),[1 4])
trig_cells=zeros(n_y-1,n_x-1,2,3);%matrix which contains trig funcions for angles in faces 2) and 4) , 
%reshape(trig_cells(5,5,:,:),[2 3])
face_centers=zeros(n_y-1,n_x-1,4,2);%matrix which contains the coordinates of points at the centres of each face
%reshape(face_centers(5,5,:,:),[4,2]) 
cell_volumes=zeros(n_y-1,n_x-1);%matrix which contains the volumen of each cell
geom_diff=zeros(n_y-1,n_x-1,4,2);%matrix for trig difussion terms(non orthogonal angles between face and nodes), not considered ghost cells 
%reshape(geom_diff(5,5,:,:),[4,2])
norm_dist_nodes=zeros(n_y-1,n_x-1,4,1);%matrix that contains the norms of the vector that join the adjacent nodes to the inner node
%reshape(norm_dist_nodes(5,5,:,:),[4,1])
wlsq_Op=zeros(n_y-1,n_x-1,2,4);%Matrix for weighted least squares operators 
%reshape(wlsq(i,j,:,:),[2,4]);
%Loop to obtain geometrical parameters 
for i=1:n_y-1
    for j=1:n_x-1
        x1=[X(i+1,j) X(i,j) X(i,j+1) X(i+1,j+1)];
        y1=[Y(i+1,j) Y(i,j) Y(i,j+1) Y(i+1,j+1)];
        %[x_cent(i,j),y_cent(i,j),len_faces(i,j,:),trig_cells(i,j,:,:),face_centers(i,j,:,:),cell_volumes(i,j),cell_centroids(i,j,:)]=cell_collocated_node(x1,y1);
        [len_faces(i,j,:),trig_cells(i,j,:,:),face_centers(i,j,:,:),cell_volumes(i,j),cell_centroids(i,j,:)]=cell_collocated_node(x1,y1);
    end
end
%___Loop to obtain geometrical parameters needed for difussion equation and
%for weighted least squares gradient 
%not included ghost cells 
for i=2:n_y-2
    for j=2:n_x-2
        nc_w=cell_centroids(i,j-1,:);%i,j-1, West
        nc_n=cell_centroids(i-1,j,:);%i-1,j, North
        nc_e=cell_centroids(i,j+1,:);%i,j+1, East
        nc_s=cell_centroids(i+1,j,:);%i+1,j, South
        nc_p=cell_centroids(i,j,:);%i,j, P
        node_cord=reshape([nc_w;nc_n;nc_e;nc_s;nc_p],[5,2]);
        angles_fns=reshape(trig_cells(i,j,:,:),[2,3]);
        adj_nodes=node_cord(1:4,:);%
        p_node=reshape(nc_p,[1,2]);%i,j
        geom_diff(i,j,:,:)= trig_diff_val(node_cord,angles_fns); 
        [wlsq_Op(i,j,:,:),norm_dist_nodes(i,j,:,:)]=gradient_lsq_weights(adj_nodes,p_node);
    end
end

%________________________Constants________________________
vel_ini=9e-5;  %29e-4; %Velocity at the inlet
rho =1.2; %Density (Kg/m3)
nu =0.0000174; %Viscosity (Kg/(m*s))
reynolds = vel_ini*rho*l_y/nu;%Reynolds number 
flux_mass=l_y*vel_ini;%mass flux
fprintf("Reynolds Number ")
disp(reynolds)

%_______________Velocity and presure fields_______________
u_vel=zeros(n_y-1,n_x-1); %Velocity in X axis
v_vel=zeros(n_y-1,n_x-1); %Velocity in Y axis
p_press=ones(n_y-1,n_x-1);  %Pressure

%________________Pressure correction coefitients___________-
d_e = zeros(size(u_vel)); % Pressure correction coefficients for the velocity field in X
d_n = zeros(size(v_vel)); % Pressure correction coefficients for the velocity field in Y

%_______________underrelaxation_factors____________________
alpha_p=0.5;
alpha_u=0.5;
alpha_v=0.5;

%________________________err________________________________
b_p =zeros(n_y-1,n_x-1);%mass balance in each cell
err=1;%error as residual of continuity
err_req=1e6; %error required
err_vec=[];%vector to storage the errors
iter=0;%iterations


%################### Main iteration #########################
%while err > err_req
%    if err > 1e12
%        fprintf("Unstable solution")
%        break
%    end
    %Momentum equations and pressure corrections
    for i=2:n_y-2
        for j =2:n_x-2
            nc_w=cell_centroids(i,j-1,:);%west node
            nc_n=cell_centroids(i-1,j,:);%north node
            nc_e=cell_centroids(i,j+1,:);%east node
            nc_s=cell_centroids(i+1,j,:);%south node 
            nc_p=cell_centroids(i,j,:);%propietary node 
            
            %adyacend node values for u
            u_w=u_vel(i,j-1);
            u_n=u_vel(i-1,j);
            u_e=u_vel(i,j+1);
            u_s=u_vel(i+1,j);
            u_p=u_vel(i,j);
            u_vec=[u_w,u_n,u_e,u_s,u_p];

            %adyacent node values for v
            v_w=u_vel(i,j-1);
            v_n=u_vel(i-1,j);
            v_e=u_vel(i,j+1);
            v_s=u_vel(i+1,j);
            v_p=u_vel(i,j);
            v_vec=[v_w,v_n,v_e,v_s,v_p];

            %aditional velocities to calculate_difussion terms
            u_nw=u_vel(i-1,j-1);
            u_ne=u_vel(i-1,j+1);
            u_sw=u_vel(i+1,j-1);
            u_se=u_vel(i+1,j+1);
            
            u_ad =[u_nw,u_ne,u_sw,u_se];
            
            v_nw=v_vel(i-1,j-1); 
            v_ne=v_vel(i-1,j+1); 
            v_sw=v_vel(i+1,j-1); 
            v_se=v_vel(i+1,j+1); 
            v_ad=[v_nw,v_ne,v_sw,v_se];

            %pressure node adjacent values 
            p_w=p_press(i,j-1);
            p_n=p_press(i-1,j);
            p_e=p_press(i,j+1);
            p_s=p_press(i+1,j);
            p_p=p_press(i,j);
            p_vec=[p_w,p_n,p_e,p_s,p_p];

            len_f = reshape(len_faces(i,j,:),[1 4]);
            angles_fns =reshape(trig_cells(i,j,:,:),[2 3]);
            diff_trig=reshape(geom_diff(i,j,:,:),[4,2]);
            delt_v=cell_volumes(i,j);
            wl_op=reshape(wlsq_Op(i,j,:,:),[2,4]);
            dist_nodes =reshape(norm_dist_nodes(i,j,:,:),[4,1]);
            [u_vel(i,j),d_e(i,j),dpy] = x_Momentum_eq(nu,rho,u_vec,v_vec,u_ad,len_f,angles_fns,diff_trig,delt_v,wl_op,p_vec,dist_nodes);
            %xmomentum
        end
    end
%{
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

