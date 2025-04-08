%function [v_star,d_Ij,dpy] = y_Momentum_eq(nu,rho,u_vec,v_vec,u_ad,len_f,angles_fns,trig_dif,delt_v,wl_op,p_vec,dist_nodes)

%adyacend node values for u
u_w=u_vec(1);
u_n=u_vec(2);
u_e=u_vec(3);
u_s=u_vec(4);
u_p=u_vec(5);

%adyacent node values for v
v_w=v_vec(1);
v_n=v_vec(2);
v_e=v_vec(3);
v_s=v_vec(4);
v_p=v_vec(5);

%aditional velocities to calculate_difussion terms
u_nw=u_ad(1);
u_ne=u_ad(2);
u_sw=u_ad(3);
u_se=u_ad(4);

% pressure at nodes W N E S and P
p_w=p_vec(1);
p_n=p_vec(2);
p_e=p_vec(3);
p_s=p_vec(4);
p_p=p_vec(5);

%lenght of faces w n, e ,s 1,2,3,4
fa_w=len_f(1);
fa_n=len_f(2);
fa_e=len_f(3);
fa_s=len_f(4);

%adyacend node values for u
u_w=u_vec(1);
u_n=u_vec(2);
u_e=u_vec(3);
u_s=u_vec(4);
u_p=u_vec(5);

%adyacent node values for v
%v_w=v_vec(1);
v_n=v_vec(2);
%v_e=v_vec(3);
v_s=v_vec(4);
v_p=v_vec(5);

%aditional velocities to calculate_difussion terms
u_nw=u_ad(1);
u_ne=u_ad(2);
u_sw=u_ad(3);
u_se=u_ad(4);

% pressure at nodes W N E S and P
p_w=p_vec(1);
p_n=p_vec(2);
p_e=p_vec(3);
p_s=p_vec(4);
p_p=p_vec(5);

%lenght of faces w n, e ,s 1,2,3,4
fa_w=len_f(1);
fa_n=len_f(2);
fa_e=len_f(3);
fa_s=len_f(4);

%adyacend node values for u
u_w=u_vec(1);
u_n=u_vec(2);
u_e=u_vec(3);
u_s=u_vec(4);
u_p=u_vec(5);

%adyacent node values for v
%v_w=v_vec(1);
v_n=v_vec(2);
%v_e=v_vec(3);
v_s=v_vec(4);
v_p=v_vec(5);

%aditional velocities to calculate_difussion terms
u_nw=u_ad(1);
u_ne=u_ad(2);
u_sw=u_ad(3);
u_se=u_ad(4);

% pressure at nodes W N E S and P
p_w=p_vec(1);
p_n=p_vec(2);
p_e=p_vec(3);
p_s=p_vec(4);
p_p=p_vec(5);

%lenght of faces w n, e ,s 1,2,3,4
fa_w=len_f(1);
fa_n=len_f(2);
fa_e=len_f(3);
fa_s=len_f(4);
%trig  for difusion terms non orthogonal angles in faces,angle phi
%trig_dif=[tan_w cos_w;tan_n cosn;tan_e cos_e;tan_s cos_s]
tanp_w=trig_dif(1,1);
cosp_w=trig_dif(1,2);
tanp_n=trig_dif(2,1);
cosp_n=trig_dif(2,2);
tanp_e=trig_dif(3,1);
cosp_e=trig_dif(3,2);
tanp_s=trig_dif(4,1);
cosp_s=trig_dif(4,2);

%distances to property node to adjacent nodes 
delta_ep_w=dist_nodes(1);
delta_ep_n=dist_nodes(2);
delta_ep_e=dist_nodes(3);
delta_ep_s=dist_nodes(4);

%delt_v volulme of the cell
%weighted least squares operator
%wl_op= [0 0 0 0;0 0 0 0];

dx=len_f(3);


Ap = []; %Coeficiente que multiplica a Up, lado izquierdo de la igualda, Salida
Anb_vab =[]; %Producto de los coeficientes Anb y velocidades Ap
S_dc =[]; %Terminos de difusión no ortogonal

%----------Difussion terms---------
%-----Face W -----
%distance between nodes  
%non ortogonal diffusion term
s_cd_w= nu*tanp_w*0.25*(v_nw + v_n - v_sw - v_s);
%diffusion term Gradient P_Ai
d_w = nu*fa_w/(cosp_w*delta_ep_w);

Ap(1)=d_w;
Anb_vab(1)=d_w*v_w;
S_dc(1)=s_cd_w;

%----- Face N -----
%non ortogonal diffusion term
s_cd_n= nu*tanp_n*0.25*(v_ne+ v_e -v_nw- v_w);
%diffusion term Gradient P_Ai
d_n = nu*fa_n/(cosp_n*delta_ep_n);

Ap(2)=d_n;
Anb_vab(2)=d_n*v_n;
S_dc(1)=s_cd_n;

%----- Face E -----
%non ortogonal diffusion term
s_cd_e= nu*tanp_e*0.25*(v_se+ v_s -v_ne- v_n);
%diffusion term Gradient P_Ai
d_e = nu*fa_e/(cosp_e*delta_ep_e);

Ap(3)=d_e;
Anb_vab(3)=d_e*v_e;
S_dc(3)=s_cd_e;

%----- Face S -----
%non ortogonal diffusion term
s_cd_s= nu*tanp_s*0.25*(v_sw+ v_w -v_se- v_e);
%diffusion term Gradient P_Ai
d_s = nu*fa_s/(cosp_s*delta_ep_s);

Ap(4)=d_s;
Anb_vab(4)=d_s*v_s;
S_dc(4)=s_cd_s;

%--------convective terms----------
%mass fluxes for U vel
%face w
fu_w=0.5*(u_w+u_p)*fa_w;
[Ap,Anb_vab]=esq_interp_upwind(fu_w,u_w,Ap,Anb_vab,fu_w);
%face n
fu_n =0.5*(u_n + u_p)*tan_fn*fa_n*rho;
%face e
fu_e = -0.5*(u_p + u_e)*fa_e*rho;
[Ap,Anb_vab]=esq_interp_upwind(fu_e,u_e,Ap,Anb_vab,fu_e);
%Face S
fu_s =-0.5*(u_p + u_s)*fa_s*tan_fs;

%mass fluxes for v
%face n
fv_n = -0.5*(v_n+v_p)*fa_n*cos_fn*rho;
[Ap,Anb_vab]=esq_interp_upwind(fv_n,u_n,Ap,Anb_vab,fv_n);
[Ap,Anb_vab]=esq_interp_upwind(fv_n,u_n,Ap,Anb_vab,fu_n);

%face S
fv_s= 0.5*(v_s+v_p)*fa_s*cos_fs*rho;
[Ap,Anb_vab]=esq_interp_upwind(fv_s,u_s,Ap,Anb_vab,fv_s);
[Ap,Anb_vab]=esq_interp_upwind(fv_s,u_s,Ap,Anb_vab,fu_s);

%Derivate of pressure X
%create vector of diference of transported quantity
delta_P = [p_w;p_n;p_e;p_s] - p_p;
grad_P= wl_op*delta_P;
dpx=grad_P(1);
dpy=grad_P(2);

%--------    Final Equation    ----------
Ap = Ap(~isnan(Ap));
Anb_vab = Anb_vab(~isnan(Anb_vab));
S_dc = S_dc(~isnan(S_dc));

a_P_i = sum(Ap); % Coefficient needed to calculate U_p and pressure correction;
d_Ij = -delt_v/ (a_P_i * dx);
u_star = (1 / a_P_i) * (sum(Anb_vab) + sum(S_dc) - dpx * delt_v); % New U_P, main output

end