function [nu_tilde_k_plus_1] = sa_transport_nutilde(nu,nu_tilde_k_vec,d_wall,delta_v,u_vec,v_vec,wl_op)
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


%____________________________________________________________

%##########       Convection     ##############  

%##########        Difusion     ############### 

%##########  Non Linear Difusion  #############
%Grad(nu_tilde_k)
delta_nu_tilde=[nu_tilde_k_w;nu_tilde_k_n;nu_tilde_k_e;nu_tilde_k_s]-nu_tilde_k_p;
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
fv2=1 - x_n/(1+x_nu*fv1);
fv2_prime=-x_nu_prime/(1+ x_nu*fv1) + (x_nu/(1 + x_nu*fv1)^2)*...
    (x_nu_prime*fv1 + x*fv1_prime);

%modified vorticity (gamma tilde)

%-------  omega_low=sqrt(omega_up_ij*omega_up_ij ---------

omega_up_12=du_dy -dv_dx;
s_vort=sqrt(2*omega_up_12^2);% vorticity (Frobenius norm of sigma_ij)
%----------------------------------------------------
s_vort_tilde= s_vort +  nu_tilde_k/((kappa^2)*(d_wall^2));%modified vorticity
s_vort_tilde_prime=(nu_tilde_k/((kappa^2)*(d_wall^2)))*fv2_prime + fv2/((kappa^2)*(d_wall^2));

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

%_Rearangement of terms 
nu_tilde_k_plus_1=(1/(delta_v*(-q_s_prime)))*(delta_v*(q_s -q_s_prime*nu_tilde_k));
end