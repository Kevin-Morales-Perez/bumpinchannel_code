function [nu_tilde_k_plus_1] = sa_transport_nutilde(nu,nu_tilde_k_vec,...
    d_wall,delta_v,u_vec,v_vec,wl_op,nu_tilde_ad,geom_disn,dist_nodes,...
    len_f,angles_fns)
%___________  Spalart Allmaras turbulence transport equation _____ 

%_____________________Constants_______________________________
kappa = 0.41; % Von Karman 
cb1=0.1355; % Basic parameter 1
cb2=0.622; %Basic parameter 2
sigma =2/3;%Prandlt Number 
cw1=cb1/kappa^2 + (1 + cb2)/sigma; %Wall parameter 1
cw2=0.3;%wall parameter 2
cw3=2;%wall parameter 3
cw3e_6=cw3^-6;%wall parameter 3 ^-6
cv1=7.1; %Viscous parameter 1 

%______      Velocity adjacent cell node values       _________
u_w=u_vec(1);
u_n=u_vec(2);
u_e=u_vec(3);
u_s=u_vec(4);
u_p=u_vec(5);

v_w=v_vec(1);
v_n=v_vec(2);
v_e=v_vec(3);
v_s=v_vec(4);
v_p=v_vec(5);
%__________________      Gradients      ______________________

%grad(U)
delta_u_vel=[u_w;u_n;u_e;u_s]-u_p;
grad_u_vel=wl_op*delta_u_vel;
%du_dx=grad_u_vel(1); Not used
du_dy=grad_u_vel(2);

%grad(V)
delta_v_vel=[v_w;v_n;v_e;v_s]-v_p;
grad_v_vel=wl_op*delta_v_vel;
dv_dx=grad_v_vel(1);
%dv_dy=grad_v_vel(2); Not used

%__________   Nu_tilde_k adjacent cell node values   ________
nu_tilde_k_w=nu_tilde_k_vec(1);
nu_tilde_k_n=nu_tilde_k_vec(2);
nu_tilde_k_e=nu_tilde_k_vec(3);
nu_tilde_k_s=nu_tilde_k_vec(4);
nu_tilde_k_p=nu_tilde_k_vec(5);

nu_tilde_k = nu_tilde_k_p;

%Additional nu_Tilde values for non orthogonal diffusion terms
nu_tilde_k_nw=nu_tilde_ad(1);
nu_tilde_k_ne=nu_tilde_ad(2);
nu_tilde_k_sw=nu_tilde_ad(3);
nu_tilde_k_se=nu_tilde_ad(4);

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

%lenght of faces w n, e ,s 1,2,3,4
fa_w=len_f(1);
fa_n=len_f(2);
fa_e=len_f(3);
fa_s=len_f(4);

%%trigonometric functions for faces 1) and 2)
sin_face_n=angles_fns(1,1);
cos_face_n=angles_fns(1,2);
%tan_face_n=angles_fns(1,3);
sin_face_s=angles_fns(2,1);
cos_face_s=angles_fns(2,2);
%tan_face_s=angles_fns(2,3);

%____________________________________________________________


%##########       Convection     ##############

Ap_c =[0 0 0 0]; %Left side contribution convective coeffitients 
Anb_nu_ab_c =[0 0 0 0]; %right side contribution convective coeffitients

%-Computing volumetric fluxes in X
fv_xw=-0.5*(u_w+u_p)*fa_w; %face W
fv_xn=-0.5*(u_n + u_p)*fa_n*sin_face_n; %face N
fv_xe=0.5*(u_p+u_e)*fa_e;% face E
fv_xs=0.5*(u_p + u_s)*sin_face_s;%face s

%-Computing volumetric fluxes in Y
%fv_yw = 0 %face W
fv_yn= 0.5*(v_n+v_p)*fa_n*cos_face_n; %face N
%fv_ye= 0 %face E
fv_ys= -0.5*(v_p + v_s)*cos_face_s;%face s

%total fluxes in face N and S
fv_n = fv_xn + fv_yn;
fv_s = fv_xs + fv_ys;

mf_omi=[fv_xw,fv_yn,fv_xe,fv_ys];%mass fluxes of most influence
phi_A_vec=[nu_tilde_k_w,nu_tilde_k_n,nu_tilde_k_e,nu_tilde_k_s];
mf_t=[fv_xw,fv_n,fv_xe,fv_s];% total mass fluxes

for i=1:4
   [Ap_c,Anb_nu_ab_c]=upwind_1st_sa(mf_omi(i),phi_A_vec(i),Ap_c,...
       Anb_nu_ab_c,mf_t(i),i);
end

%##########        Difusion     ############### 

Ap_d = []; %LSE diffusion
Anb_uab_d =[]; %RSE  diffusion 
S_dc =[]; %non orthogonal diffusion terms
gamma_up_sa=nu_tilde_k + nu;%coefficient for linear diffusion in SA eq.

%-----Face W -----
%non ortogonal diffusion term
s_cd_w= gamma_up_sa*tanp_w*0.25*(nu_tilde_k_nw + nu_tilde_k_n - ...
    nu_tilde_k_sw - nu_tilde_k_s);
%diffusion term Gradient P_Ai
d_w = gamma_up_sa*fa_w/(cosp_w*delta_ep_w);

Ap_d=[Ap_d,d_w];
Anb_uab_d= [Anb_uab_d,d_w*u_w];
S_dc(1)=s_cd_w;

%----- Face N -----
%non ortogonal diffusion term
s_cd_n= gamma_up_sa*tanp_n*0.25*(nu_tilde_k_ne+ nu_tilde_k_e -...
    nu_tilde_k_nw- nu_tilde_k_w);
%diffusion term Gradient P_Ai
d_n = gamma_up_sa*fa_n/(cosp_n*delta_ep_n);

Ap_d=[Ap_d,d_n];
Anb_uab_d=[Anb_uab_d,d_n*u_n];
S_dc(2)=s_cd_n;

%----- Face E -----
%non ortogonal diffusion term
s_cd_e= gamma_up_sa*tanp_e*0.25*(nu_tilde_k_se+ nu_tilde_k_s-...
    nu_tilde_k_ne- nu_tilde_k_n);
%diffusion term Gradient P_Ai
d_e = gamma_up_sa*fa_e/(cosp_e*delta_ep_e);

Ap_d=[Ap_d,d_e];
Anb_uab_d=[Anb_uab_d,d_e*u_e];
S_dc(3)=s_cd_e;

%----- Face S -----
%non ortogonal diffusion term
s_cd_s= gamma_up_sa*tanp_s*0.25*(nu_tilde_k_sw+ nu_tilde_k_w...
    -nu_tilde_k_se- nu_tilde_k_e);
%diffusion term Gradient P_Ai
d_s = gamma_up_sa*fa_s/(cosp_s*delta_ep_s);

Ap_d=[Ap_d,d_s];
Anb_uab_d=[Anb_uab_d,d_s*u_s];
S_dc(4)=s_cd_s;


%##########  Non Linear Difusion  #############
%Grad(nu_tilde_k)
delta_nu_tilde=[nu_tilde_k_w;nu_tilde_k_n;nu_tilde_k_e;nu_tilde_k_s]...
    -nu_tilde_k_p;
grad_nu_tilde=wl_op*delta_nu_tilde;
dnu_tilde_dx=grad_nu_tilde(1);
dnu_tilde_dy=grad_nu_tilde(2);

nonlin_diff= (cb2/sigma)*(dnu_tilde_dx^2 + dnu_tilde_dy^2);

%##########        Sources        #############

%  Implicit Taylor 1st order expansion nu_tilde_k+1

%Q(nu_tilde)=Tp(Nu_tilde) + Td(nu_tilde)

%Q(nu_tilde_k+1)=Q(nu_tilde_k)+Q'(nu_tilde_k)*(nu_tilde_k+1-nu_tilde_k)

%Q'(nu_tilde_k)=Tp'(nu_tilde_k) + Td'(nu_tilde_k)

%__________ Turbulence Production Tp(nu_tilde_k) ___________

%viscous ratio x_nu
x_nu=nu_tilde_k/nu;
x_nu_prime=1/nu;

%Viscous damping function 1 fv1
fv1= x_nu^3/(x_nu^3 +cv1^3);
fv1_prime = ((3*x_nu_prime^2)/(x_nu^3+ cv1^3))*x_nu_prime - x_nu_prime*...
    ((3*(x_nu^5))/((x_nu^3 +cv1^3)^2));

%Viscous damping function 2 fv2
fv2=1 - x_nu/(1+x_nu*fv1);
fv2_prime=-x_nu_prime/(1+ x_nu*fv1) + (x_nu/(1 + x_nu*fv1)^2)*...
    (x_nu_prime*fv1 + x_nu*fv1_prime);

%modified vorticity (gamma tilde)

%-------  omega_low=sqrt(omega_up_ij*omega_up_ij ---------

omega_up_12=du_dy -dv_dx;
s_vort=sqrt(2*omega_up_12^2);% vorticity (Frobenius norm of sigma_ij)
%----------------------------------------------------
s_vort_tilde= s_vort +  nu_tilde_k/((kappa^2)*(d_wall^2));%modified vorticity
s_vort_tilde_prime=(nu_tilde_k/((kappa^2)*(d_wall^2)))*fv2_prime +...
    fv2/((kappa^2)*(d_wall^2));

%Tp
tp=cb1*s_vort_tilde*nu_tilde_k;
tp_prime=cb1*s_vort_tilde + cb1*s_vort_tilde_prime*nu_tilde_k;

%__________ Turbulence Destruction Td(nu_tilde_k) __________

%function r
r=nu_tilde_k/(s_vort_tilde*(kappa^2)*(d_wall^2));
r_prime=1/(s_vort_tilde*(kappa^2)*(d_wall^2))-(nu_tilde_k/(s_vort*...
    (kappa^2)*(d_wall^2) + nu_tilde_k*fv2))*(fv2 + nu_tilde_k*fv2_prime);

%function g
g=r + cw2*(r^6 -r);
g_prime=(1+cw2*(6*r^5 - 1))*r_prime;

%wall damping function fw
fw=g*((1+cw3^6)/(g^6 + cw3^6))^(1/6);
fw_prime=(-1/6)*(((cw3e_6 + g^-6)/(cw3e_6 + 1))^(-7/6))*((-6*g^-7)/...
    (1+cw3e_6))*g_prime;

%Td
td=-cw1*fw*(nu_tilde_k/d_wall)^2;
td_prime=-cw1*fw_prime*(nu_tilde_k/d_wall)^2 - cw1*fw*2*nu_tilde_k/...
    (d_wall^2);

q_s=tp + td;
q_s_prime=tp_prime + td_prime;

%###################### Final Subtitution #############################
%--------    Final Equation    ----------
Ap_c = Ap_c(~isnan(Ap_c));
Anb_nu_ab_c= Anb_nu_ab_c(~isnan(Anb_nu_ab_c));

Ap_d = Ap_d(~isnan(Ap_d));
Anb_uab_d = Anb_uab_d(~isnan(Anb_uab_d));
S_dc = S_dc(~isnan(S_dc));

Ap_c=sum(Ap_c);
Anb_nu_ab_c=sum(Anb_nu_ab_c);

Ap_d=sum(Ap_d);
Anb_uab_d=sum(Anb_uab_d);
S_dc=sum(S_dc);


left_s_coeff= Ap_c + Ap_d - delta_v*q_s_prime;
right_s_coeff= Anb_nu_ab_c + Anb_uab_d + S_dc+ delta_v*(nonlin_diff + q_s - q_s_prime*nu_tilde_k);

%_Rearangement of terms 
nu_tilde_k_plus_1=(1/left_s_coeff)*right_s_coeff;
end