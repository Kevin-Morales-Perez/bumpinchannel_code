function [u_star,d_Ij,dpy,dx_tau_xy] = x_Momentum_eq(nu,rho,u_vec,v_vec,u_ad,len_f,angles_fns,geom_disn,delt_v,wl_op,p_vec,tau_xx_vec,tau_xy_vec,dist_nodes)
% x Momentum equation 
% Finite volume x Momentum equation
%Example of imput values
%nu=1; %viscosity
%rho=1;% density
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

%tau_xx at nodes W N E S and P
tau_xx_w=tau_xx_vec(1);
tau_xx_n=tau_xx_vec(2);
tau_xx_e=tau_xx_vec(3);
tau_xx_s=tau_xx_vec(4);
tau_xx_p=tau_xx_vec(5);

%tau_xy at nodes W N E S and P
tau_xy_w=tau_xy_vec(1);
tau_xy_n=tau_xy_vec(2);
tau_xy_e=tau_xy_vec(3);
tau_xy_s=tau_xy_vec(4);
tau_xy_p=tau_xy_vec(5);

%lenght of faces w n, e ,s 1,2,3,4
fa_w=len_f(1);
fa_n=len_f(2);
fa_e=len_f(3);
fa_s=len_f(4);

%trigonometric functions for faces 1) and 2)
tan_fn=angles_fns(1,3);
tan_fs=angles_fns(2,3);
cos_fn=angles_fns(1,2);
cos_fs=angles_fns(2,2);

%trig  for difusion terms non orthogonal angles in faces,angle phi
%geom_disn=[tan_w cos_w;tan_n cosn;tan_e cos_e;tan_s cos_s]
tanp_w =geom_disn(1,1);
cosp_w =geom_disn(1,2);
tanp_n =geom_disn(2,1);
cosp_n =geom_disn(2,2);
tanp_e =geom_disn(3,1);
cosp_e =geom_disn(3,2);
tanp_s =geom_disn(4,1);
cosp_s =geom_disn(4,2);

%distances to property node to adjacent nodes 
delta_ep_w=dist_nodes(1);
delta_ep_n=dist_nodes(2);
delta_ep_e=dist_nodes(3);
delta_ep_s=dist_nodes(4);

%delt_v volulme of the cell
%weighted least squares operator
%wl_op= [0 0 0 0;0 0 0 0];

%_________________________________________________

dx=len_f(4);

%----------coeficientes Ap -  Anb - Scd ---------

Ap = []; %LSE
Anb_uab =[]; %RSE 
S_dc =[]; %non orthogonal diffusion terms

%----------Difussion terms---------
%-----Face W -----
%non ortogonal diffusion term
s_cd_w= nu*tanp_w*0.25*(u_nw + u_n - u_sw - u_s);
%diffusion term Gradient P_Ai
d_w = nu*fa_w/(cosp_w*delta_ep_w);

Ap(1)=d_w;
Anb_uab(1)=d_w*u_w;
S_dc(1)=s_cd_w;

%----- Face N -----
%non ortogonal diffusion term
s_cd_n= nu*tanp_n*0.25*(u_ne+ u_e -u_nw- u_w);
%diffusion term Gradient P_Ai
d_n = nu*fa_n/(cosp_n*delta_ep_n);

Ap(2)=d_n;
Anb_uab(2)=d_n*u_n;
S_dc(1)=s_cd_n;

%----- Face E -----
%non ortogonal diffusion term
s_cd_e= nu*tanp_e*0.25*(u_se+ u_s -u_ne- u_n);
%diffusion term Gradient P_Ai
d_e = nu*fa_e/(cosp_e*delta_ep_e);

Ap(3)=d_e;
Anb_uab(3)=d_e*u_e;
S_dc(3)=s_cd_e;

%----- Face S -----
%non ortogonal diffusion term
s_cd_s= nu*tanp_s*0.25*(u_sw+ u_w -u_se- u_e);
%diffusion term Gradient P_Ai
d_s = nu*fa_s/(cosp_s*delta_ep_s);

Ap(4)=d_s;
Anb_uab(4)=d_s*u_s;
S_dc(4)=s_cd_s;

%--------convective terms----------
%mass fluxes for U vel
%face w
fu_w=0.5*(u_w+u_p)*fa_w*rho;
[Ap,Anb_uab]=esq_interp_upwind(fu_w,u_w,Ap,Anb_uab,fu_w);
%face n
fu_n =0.5*(u_n + u_p)*tan_fn*fa_n*rho;
%face e
fu_e = -0.5*(u_p + u_e)*fa_e*rho;
[Ap,Anb_uab]=esq_interp_upwind(fu_e,u_e,Ap,Anb_uab,fu_e);
%Face s
fu_s =-0.5*(u_p + u_s)*fa_s*tan_fs*rho;

%mass fluxes for v
%face n
fv_n = -0.5*(v_n+v_p)*fa_n*cos_fn*rho;
[Ap,Anb_uab]=esq_interp_upwind(fv_n,u_n,Ap,Anb_uab,fv_n);
[Ap,Anb_uab]=esq_interp_upwind(fv_n,u_n,Ap,Anb_uab,fu_n);

%face S
fv_s= 0.5*(v_s+v_p)*fa_s*cos_fs*rho;
[Ap,Anb_uab]=esq_interp_upwind(fv_s,u_s,Ap,Anb_uab,fv_s);
[Ap,Anb_uab]=esq_interp_upwind(fv_s,u_s,Ap,Anb_uab,fu_s);

%Derivate of pressure X
%create vector of diference of transported quantity
delta_P = [p_w;p_n;p_e;p_s] - p_p;
grad_P= wl_op*delta_P;
dpx=grad_P(1);
dpy=grad_P(2);


%Derivative of u'2- respect X
delta_tau_xx=[tau_xx_w;tau_xx_n;tau_xx_e;tau_xx_s]-tau_xx_p;
grad_tau_xx=wl_op*delta_tau_xx;
dx_tau_xx=grad_tau_xx(1);
%dy_tau_xx=grad_tau_xx(2);

%Derivative of uv'- respect y
delta_tau_xy=[tau_xy_w;tau_xy_n;tau_xy_e;tau_xy_s]-tau_xy_p;
grad_tau_xy=wl_op*delta_tau_xy;
dx_tau_xy=grad_tau_xy(1);
dy_tau_xy=grad_tau_xy(2);

%--------    Final Equation    ----------
Ap = Ap(~isnan(Ap));
Anb_uab = Anb_uab(~isnan(Anb_uab));
S_dc = S_dc(~isnan(S_dc));

a_P_i = sum(Ap); % Coefficient needed to calculate U_p and pressure correction;
d_Ij = -delt_v/ (a_P_i * dx);
u_star = (1 / a_P_i) * (sum(Anb_uab) + sum(S_dc) - dpx * delt_v - dx_tau_xx*delt_v - dy_tau_xy*delt_v); % New U_P, main output

end