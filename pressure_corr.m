function [p_star] = pressure_corr(u_vec,v_vec,p_vec,d_e,d_n)
%Pressure correction
%Pressure correction for collocated grid for SIMPLE alghoritm


%nodal u velocities 
u_w=u_vec(1);
u_n=u_vec(2);
u_e=u_vec(3);
u_s=u_vec(4);
u_p=u_vec(5);
u_ee=u_vec(6);

%nodal v velocities 
v_w=v_vec(1);
v_n=v_vec(2);
v_e=v_vec(3);
v_s=v_vec(4);
v_p=v_vec(5);
v_ee=v_vec(6);

%nodal presures 
p_w=p_vec(1);
p_n=p_vec(2);
p_e=p_vec(3);
p_s=p_vec(4);
p_p=p_vec(5);
p_ee=p_vec(6);


u_ef= 0.5*(up+ue) + 0.5*(dp + de)*(pp-pe) - 0.25*dp*(pw-pe) - 0.25*de*(pp-pee);
%compute velocities at the faces using Rie Chow formula
%compute mass fluxes
%use the previos mass fluxes to compute the discretized continuity equation
%for pressure correction
%get pressure corrected
end