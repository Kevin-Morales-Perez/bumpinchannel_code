function [nu_tilde_k_plus_1] = sa_transport_nutilde(nu,nu_tilde_k,d_wall,delta_v)
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

%##########       Convection     ##############  

%##########        Difusion     ############### 

%##########  Non Linear Difusion  #############

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
fv1= x_nu^3/(x_nu^3 +cv1);
fv1_prime = ((3*x_nu_prime^2)/(x_nu^3+ cv1^3))*x_nu_prime - x_nu_prime*...
    ((3*x_nu^5)/(x_nu^3 +cv1^3));

%Viscous damping function 2 fv2
fv2=1 - x_n/(1+x_nu*fv1);
fv2_prime=-x_nu_prime/(1+ x_nu*fv1) + (x_nu/(1 + x_nu*fv1)^2)*(x_nu_prime*fv1 + x*fv1_prime);

%modified vorticity (gamma tilde)
gamma=0;%Frobenius norm of rotation tensor or vorticity (Pending)
gamma_tilde= gamma +  nu_tilde_k/((kappa^2)*(d_wall^2));%modified vorticity
gamma_tilde_prime=(nu_tilde_k/((kappa^2)*(d_wall^2)))*fv2_prime + fv2/((kappa^2)*(d_wall^2));

%Tp
tp=cb1*gamma_tilde*nu_tilde_k;
tp_prime=cb1*gamma_tilde + cb1*gamma_tilde_prime*nu_tilde_k;

%__________ Turbulence Destruction Td(nu_tilde_k) __________

%function r
r=nu_tilde_k/(gamma_tilde*kappa^2*d_wall^2);
r_prime=1/(gamma_tilde*kappa^2*d_wall^2) - (nu_tilde_k/(gamma*kappa^2*...
    d_wall^2 + nu_tilde_k*fv2))*(fv2 + nu_tilde_k*fv2_prime);

%function g
g=r + cw2*(r^6 -r);
g_prime=(1+cw2*(6*r^5 - 1))*r_prime;

%wall damping function fw
fw=g*((1+cw3^6)/(g^6 + cw3^6))^(1/6);
fw_prime=(-1/6)*(((cw3e_6 + g^-6)/(cw3e_6 + 1))^(7/6))*((-6*g^-7)/(1+cw3e_6*g_prime;

%Td
td=-cw1*fw*(nu_tilde_k/d_wall)^2;
td_prime=-cw1*fw_prime*(nu_tilde_k/d_wall)^2 - cw1*fw*2*nu_tilde_k/(d_wall^2);

q_s=tp + td;
q_s_prime=tp_prime + td_prime;

nu_tilde_k_plus_1=(1/(delta_v*(-q_s_prime)))*(delta_v*(q_s -q_s_prime*nu_tilde_k));
end