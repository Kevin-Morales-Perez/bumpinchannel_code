function [nu_tilde_k_f] = sa_transport_nutilde(nu,rho,nu_tilde_k)
%Spalart Allmaras turbulence transpor equation 


%_____________________Constants_______________________________
kappa = 0.41; % Von Karman 
cb1=0.1355; % Basic parameter 1
cb2=0.622; %Basic parameter 2
sigma =2/3;%Prandlt Number 
cw1=cb1/kappa^2 + (1 + cb2)/sigma; %Wall parameter 1
cw2=0.3;%wall parameter 2
cw3=2;%wall parameter 3
cw3e_6=cw3^-6;
cv1=7.1; %Viscous parameter 1 

%___________      Convection      ___________          
%__________        Difusion     _____________                           
%__________  Non Linear Difusion  ___________          
%__________ Turbulence Production ___________           
%__________ Turbulence Destruction __________
%Implicit (Taylor 1st order expansion of fw)




X_diff_k=(1/nu)*(nu_tilde_diff);

fv1_diff_k=((3*X_k^2)/(X_k^3 +cv1^3))*X_diff_k - (X_k^3)/(X_k^3+cv1^3)*...
    3*(X_k^2)*X_diff_k;

fv2_diff_k=-X_diff_k/(1+ X_k*fv1_k) + X_k/(1+X_k*fv1_k)^2)*(X_k*fv1_diff_k ...
    + X_k*fv1_diff_k;

r_diff_k= (nu_tilde_diff/(gamma_tilde*(kappa^2)*(d_wall^2)))  + ... 
    (nu_tilde_k/(vorticity_m*(kappa^2)*d_wall^2 + ...
    nu_tilde_k*fv2_k)^2)*(nu_tilde_diff*fv2_k + nu_tilde_k*fv2_diff_k); 

g_diff_k=1+cw2*(r_k^6 -r_k)*r_diff_k;

f_w_diff_k=-1/6*(((cw3e_6+g_k^-6)/(cw3e_6+1))^-(7/6))*((-6*g_k^-7)/...
    (1+cw3e_6))*g_diff_k;



end