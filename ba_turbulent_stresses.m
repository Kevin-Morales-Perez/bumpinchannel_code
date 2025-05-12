function [tau_xx_val,tau_xy_val,tau_yy_val] = ba_turbulent_stresses(u_vec,v_vec,wl_op,nu_turbulent_val)
%Computation of Reynolds Stress tensor using Bousinessq asumption

%adyacent node values for u
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

%Computing gradients

%grad(U)
delta_u_vel=[u_w;u_n;u_e;u_s]-u_p;
grad_u_vel=wl_op*delta_u_vel;
du_dx=grad_u_vel(1);
du_dy=grad_u_vel(2);


%grad(V)
delta_v_vel=[v_w;v_n;v_e;v_s]-v_p;
grad_v_vel=wl_op*delta_v_vel;
dv_dx=grad_v_vel(1);
dv_dy=grad_v_vel(2);


tau_xx_val=-2*nu_turbulent_val*du_dx;
tau_xy_val=-nu_turbulent_val*(du_dy+dv_dx);
tau_yy_val=-2*nu_turbulent_val*dv_dy;


end