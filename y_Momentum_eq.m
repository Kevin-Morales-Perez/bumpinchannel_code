function [v_star,d_iJ] = y_Momentum_eq(nu,rho,u_vec,v_vec,v_ad,len_f,angles_fns,geom_disn,delt_v,wl_op,tau_yy_vec,dist_nodes,dpy,dx_tau_xy)
% y Momentum equation 
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
v_w=v_vec(1);
v_n=v_vec(2);
v_e=v_vec(3);
v_s=v_vec(4);
v_p=v_vec(5);

%aditional velocities to calculate_difussion terms
v_nw=v_ad(1);
v_ne=v_ad(2);
v_sw=v_ad(3);
v_se=v_ad(4);

%tau_yy at nodes W N E S  and P
tau_yy_w=tau_yy_vec(1);
tau_yy_n=tau_yy_vec(2);
tau_yy_e=tau_yy_vec(3);
tau_yy_s=tau_yy_vec(4);
tau_yy_p=tau_yy_vec(5);

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

dy=0.5*(fa_w+fa_e);

%----------Initializing coefficients Ap - Anb - Scd ---------

Ap = []; % Coefficient that multiplies Up, left side of the equation, Outlet
Anb_vab = []; % Product of coefficients Anb and velocities Ap
s_dc = []; % Non-orthogonal diffusion terms

%----------Diffusion Terms---------
%-------------Face W -------------
% Non-orthogonal diffusion term 
s_cd_w =nu*tanp_w*0.25*(v_nw+ v_n-v_sw- v_s);
% Diffusion terms Pressure Gradient P_Ai
d_w = nu*fa_w/(cosp_w*delta_ep_w);

Ap(1)=d_w;
Anb_vab(1)=d_w*v_w;
s_dc(1)=s_cd_w;

%----- Face N -----
% Non-orthogonal diffusion term 
s_cd_n = nu*tanp_n*0.25*(v_ne+ v_e-v_nw- v_w);
% Diffusion terms Pressure Gradient P_Ai
d_n = nu*fa_n/(cosp_n*delta_ep_n);

Ap(2)=d_n;
Anb_vab(2)=d_n*v_n;
s_dc(2)= s_cd_n ;

%----- Face E ------
% Non-orthogonal diffusion term 
s_cd_e = nu*tanp_e*0.25*(v_se+ v_s-v_ne- v_n);
% Diffusion terms Pressure Gradient P_Ai
d_e = nu*fa_e/(cosp_e*delta_ep_e);

Ap(3)=d_e;
Anb_vab(3)=d_e*v_e;
s_dc(3)= s_cd_e;

% Face S
% Non-orthogonal diffusion term 
s_cd_s= nu*tanp_s*0.25*(v_sw+ v_w-v_se- v_e);

% Diffusion terms Pressure Gradient P_Ai
d_s = nu*fa_s/(cosp_s*delta_ep_s);

Ap(4)=d_s;
Anb_vab(4)=d_s*v_s;
s_dc(4)= s_cd_s;


%--------Convective terms----------

%Mass fluxes for  U ,transported variable v
%-------------Face W -------------
Fu_w=0.5*(u_w+u_p)*fa_w*rho;
[Ap,Anb_vab]=esq_interp_upwind(Fu_w,v_w,Ap,Anb_vab,Fu_w);

%-------------Face N -------------
Fu_n =0.5*(u_n+u_p)*fa_n*cos_fn*tan_fn*rho;

%-------------Face E -------------
Fu_e = -0.5*(u_e+u_p)*fa_e*rho;
[Ap,Anb_vab]=esq_interp_upwind(Fu_e,v_e,Ap,Anb_vab,Fu_e);

%-------------Face S ------------
Fu_s = -0.5*(u_s+u_p)*fa_s*tan_fs*cos_fs*rho;


%Mass fluxes for  V ,transported variable v

%-------------Face W -------------
%Fv_w=0;

%-------------Face N -------------
Fv_n =-0.5*(v_n + v_p)*fa_n*cos_fn*rho;
[Ap,Anb_vab]=esq_interp_upwind(Fv_n,v_n,Ap,Anb_vab,Fv_n);
[Ap,Anb_vab]=esq_interp_upwind(Fv_n,v_n,Ap,Anb_vab,Fu_n);


%-------------Face E -------------
%Fv_E=0;

%-------------Face S ------------
Fv_s=0.5*(v_s+v_s)*fa_s*cos_fs*rho;
[Ap,Anb_vab]=esq_interp_upwind(Fv_s,v_s,Ap,Anb_vab,Fv_s);
[Ap,Anb_vab]=esq_interp_upwind(Fv_s,v_s,Ap,Anb_vab,Fu_s);

%Derivate of pressure respect y
%obtained from X momentum equation

%Derivative of uv'- respect X
%Obtained from  X momentum equation

%Derivative of v2'- respect y
delta_tau_yy=[tau_yy_w;tau_yy_n;tau_yy_e;tau_yy_s]-tau_yy_p;
grad_tau_yy=wl_op*delta_tau_yy;
%dx_tau_yy=grad_tau_yy(1);
dy_tau_yy=grad_tau_yy(2);


%----------------------------- ----------
%--------    Final Equation    ----------
% Final coefficient AP, left side of the equality;

Ap=Ap(~isnan(Ap));
Anb_vab=Anb_vab(~isnan(Anb_vab));
s_dc=s_dc(~isnan(s_dc));

A_P_i = sum(Ap); % Coefficient needed to calculate U_p and pressure correction
d_iJ=-delt_v/(A_P_i*dy); %coeficient pressure corrections  
v_star =(1/A_P_i)*(sum(Anb_vab) + sum(s_dc) - dpy*delt_v - dx_tau_xy*delt_v - dy_tau_yy*delt_v); %New U_P, main output

end






















