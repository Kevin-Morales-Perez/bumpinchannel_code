%Bump in channel problem
%Written by Kevin Morales 
%Instituto Politécnico Nacional , Aeronautical Engineering
%clear all
close all 
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
%reshape(cell_centroids(i,j,:),[1,2])
len_faces=zeros(n_y-1,n_x-1,4); %matrix which contains lenght of the 4 faces of each cell, 
%reshape(len_faces(i,j,:),[1 4])
uvec_norm_faces=zeros(n_y-1,n_x-1,4,2);% Matrix with unitary vectors normal to the faces
%reshape(uvec_dist_nodes(i,j,:,:),[4,2])
trig_cells=zeros(n_y-1,n_x-1,2,3);%matrix which contains trig funcions for angles in faces 2) and 4) , 
%reshape(trig_cells(i,j,:,:),[2 3])
face_centers=zeros(n_y-1,n_x-1,4,2);%matrix which contains the coordinates of points at the centres of each face
%reshape(face_centers(i,j,:,:),[4,2]) 
cell_volumes=zeros(n_y-1,n_x-1);%matrix which contains the volumen of each cell
geom_diff=zeros(n_y-1,n_x-1,4,2);%matrix for trig difussion terms(non orthogonal angles between face and nodes), not considered ghost cells 
%reshape(geom_diff(i,j,:,:),[4,2])
uvec_dist_nodes=zeros(n_y-1,n_x-1,4,2);%matrix that contains the unitary vectors that join the adjacent nodes to the inner node
%reshape(uvec_dist_nodes(i,j,:,:),[4,2])
norm_dist_nodes=zeros(n_y-1,n_x-1,4,1);%matrix that contains the norms of the vectors that join the adjacent nodes to the inner node
%reshape(norm_dist_nodes(i,j,:,:),[4,1])
wlsq_Op=zeros(n_y-1,n_x-1,2,4);%Matrix for weighted least squares operators 
%reshape(wlsq_Op(i,j,:,:),[2,4]);
d_mnw=zeros(n_y-1,n_x-1);%Minimun distance to the nearest wall for Spalart - Allmaras Equation
%Loop to obtain geometrical parameters , includet ghost cells
for i=1:n_y-1
    for j=1:n_x-1
        x1=[X(i+1,j) X(i,j) X(i,j+1) X(i+1,j+1)];
        y1=[Y(i+1,j) Y(i,j) Y(i,j+1) Y(i+1,j+1)];
         [len_faces(i,j,:),trig_cells(i,j,:,:),face_centers(i,j,:,:),cell_volumes(i,j),cell_centroids(i,j,:),uvec_norm_faces(i,j,:,:)]=cell_collocated_node(x1,y1);
    end
end
%___Loop to obtain geometrical parameters needed for difussion equation and
%for weighted least squares gradient, minimum distance to the wall
%for Spalart - Allmaras 
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
        [wlsq_Op(i,j,:,:),norm_dist_nodes(i,j,:,:),uvec_dist_nodes(i,j,:,:)]=gradient_lsq_weights(adj_nodes,p_node);
    end
end

%minimum distance to the wall for Spalart - Allmaras eq.
%not included ghost cells 
for i=2:n_y-2
    k=0;
    for j=2:(0.5*n_x)
        jinv=n_x-2-k;
        k=k+1;
        nod_cord=reshape(cell_centroids(i,j,:),[1,2]);
        d_mnw(i,j)=near_walld(nod_cord);
        d_mnw(i,jinv)=d_mnw(i,j);
    end
end
clear k;


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
p_corr=zeros(n_y-1,n_x-1); %pressure corrections

%________________Pressure correction coefitients___________-
d_e = zeros(size(u_vel)); % Pressure correction coefficients for the velocity field in X
d_n = zeros(size(v_vel)); % Pressure correction coefficients for the velocity field in Y

%_______________derivative of  pressure respect y___________
dp_dy=zeros(n_y-1,n_x-1);
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
            
            %pressure node adjacent values 
            p_w=p_press(i,j-1);
            p_n=p_press(i-1,j);
            p_e=p_press(i,j+1);
            p_s=p_press(i+1,j);
            p_p=p_press(i,j);
            p_vec=[p_w,p_n,p_e,p_s,p_p];

            len_f = reshape(len_faces(i,j,:),[1 4]);
            angles_fns =reshape(trig_cells(i,j,:,:),[2 3]);
            geom_disn=reshape(geom_diff(i,j,:,:),[4,2]);
            delt_v=cell_volumes(i,j);
            wl_op=reshape(wlsq_Op(i,j,:,:),[2,4]);
            dist_nodes =reshape(norm_dist_nodes(i,j,:,:),[4,1]);
            [u_vel(i,j),d_e(i,j),dp_dy(i,j)] = x_Momentum_eq(nu,rho,u_vec,v_vec,u_ad,len_f,angles_fns,geom_disn,delt_v,wl_op,p_vec,dist_nodes);
            %xmomentum
        end
    end
    
   %x boundary conditions
    for i=2:n_y-2
            for j =2:n_x-2
                
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
                v_nw=v_vel(i-1,j-1); 
                v_ne=v_vel(i-1,j+1); 
                v_sw=v_vel(i+1,j-1); 
                v_se=v_vel(i+1,j+1); 
                v_ad=[v_nw,v_ne,v_sw,v_se];
    
                len_f = reshape(len_faces(i,j,:),[1 4]);
                angles_fns =reshape(trig_cells(i,j,:,:),[2 3]);
                geom_disn=reshape(geom_diff(i,j,:,:),[4,2]);
                delt_v=cell_volumes(i,j);
                dist_nodes =reshape(norm_dist_nodes(i,j,:,:),[4,1]);
                dpy=dp_dy(i,j);
                [v_vel(i,j),d_n(i,j)] = y_Momentum_eq(nu,rho,u_vec,v_vec,v_ad,len_f,angles_fns,geom_disn,delt_v,dist_nodes,dpy);
                %xmomentum
            end
   end
   
    %y_boundary
    %pressure correction
    for i=3:n_y-3
            for j =3:n_x-3
                %Matrix for pressures
                p_m=[0,0,p_press(i-2,j),0,0;0,p_press(i-1,j-1),p_press(i-1,j)...
                    ,p_press(i-1,j+1),0;p_press(i,j-2),p_press(i,j-1),...
                    p_press(i,j),p_press(i,j+1),p_press(i,j+2);0,...
                    p_press(i+1,j-1),p_press(i+1,j),p_press(i+1,j+1),0;...
                    0,0,p_press(i+2,j),0,0];
                %lenght of the faces
                len_faces_cell=reshape(len_faces(i,j,:),[1 4]);
                %Distances from node P to nodes W,N,S,E
                d_eps_vec=reshape(norm_dist_nodes(i,j,:,:),[4,1])';
                %Unitary vectors from node P to nodes W,N,S,E
                eps_vec=reshape(uvec_dist_nodes(i,j,:,:),[4,2]);
                %Unitary vectors normal to faces
                normf_vec=reshape(uvec_dist_nodes(i,j,:,:),[4,2]);
                %vectors with nodal velocities
                u_W=[u_vel(i,j-1),v_vel(i,j-1)];
                u_N=[u_vel(i-1,j),v_vel(i-1,j)];
                u_E=[u_vel(i,j-1),v_vel(i,j-1)];
                u_S=[u_vel(i+1,j),v_vel(i+1,j)];
                u_P=[u_vel(i,j),v_vel(i,j)];
                u_vec=[u_W;u_N;u_E;u_S;u_P];
                %coeffitients
                d_k=[1,1,1,1];%Not yet defined
                %weighted least squares operator matrix
                wl_op_w=reshape(wlsq_Op(i,j-1,:,:),[2,4]);
                wl_op_n=reshape(wlsq_Op(i-1,j,:,:),[2,4]);
                wl_op_e=reshape(wlsq_Op(i,j+1,:,:),[2,4]);
                wl_op_s=reshape(wlsq_Op(i+1,j,:,:),[2,4]);
                wl_op_p=reshape(wlsq_Op(i,j,:,:),[2,4]);
                wl_op_mat=[wl_op_w;wl_op_n;wl_op_e;wl_op_s;wl_op_p];
                %pressure correction function
                p_corr(i,j)=pressure_corr(p_m,len_faces_cell,d_eps_vec,eps_vec,normf_vec,u_vec,d_k,wl_op_mat);
            end
    end

    %{
    %calculating err
end
   %}

toc

