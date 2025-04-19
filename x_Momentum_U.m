%function[u_star,d_Ij,dpy]=x_Momentum_eq(cent_vec,u_vec,v_vec,u_ad,v_ad,p_vec,len_f,angles_fns,diff_trig)
%base for function equation in x
%base script for x momentum equation , when finished it will be coverted to
%a function
%input  values 
nu=1; %viscosity
rho=1;% density

%nodes coordinates W, N , E, S , (p centroid of the cell)
nc_w =[0 0];
nc_n =[0 0];
nc_e =[0 0];
nc_s= [0 0];
nc_p= [0 0];

%adyacend node values for u
u_w=1;
u_n=1;
u_e=1;
u_s=1;
u_p=1;

u_vec=[u_w,u_n,u_e,u_s,u_p];

%adyacent node values for v
v_w=1;
v_n=1;
v_e=1;
v_s=1;
v_p=1;


v_vec=[u_w,u_n,u_e,u_s,u_p];

%aditional velocities to calculate_difussion terms
u_nw=0;
u_ne=0;
u_sw=0;
u_se=0;

u_ad =[u_nw,u_ne,u_sw,u_se];

v_nw=0; 
v_ne=0; 
v_sw=0; 
v_se=0; 
v_ad=[v_nw,v_ne,v_sw,v_se];

%lenght of faces w n, e ,s 1,2,3,4
len_f=[0 0 0 0];

%trigonometric funtions for faces 1 ,N) and 2 ,S)
angles_fns= [0 0 0;0 0 0];

%angle phi tangents ( used for difussion eq)
tan_w=(nc_w(2)-nc_p(2))/(nc_w(1)-nc_p(1));
phi_w=atan(tan_w);
cos_w=cos(phi_w);

tan_n=tan(atan(angles_fns(1,1)) + atan2(nc_n(1)-nc_p(1),nc_n(2)-nc_p(2)));
phi_n=atan(tan_n);
cos_n=cos(phi_n);

tan_e=(nc_p(2)-nc_p(2))/(nc_e(1)-nc_p(1));
phi_e=atan(tan_e);
cos_e=cos(phi_e);

tan_s =tan(atan(angles_fns(2,1)) +  atan2(nc_p(1)-nc_s(1),nc_p(2)-nc_s(2)));
phi_s=atan(tan_s);
cos_s=cos(phi_s);

%trig funtions for non orthogonal terms


delt_v=1;%volume of the cell

%weighted least squares operator
wl_op= [0 0 0 0;0 0 0 0];

%Get the values of pressure at nodes W N E S and P
p_w=0;
p_n=0;
p_e=0;
p_s=0;
p_p=0;
%_________________________________________________

dx=len_f(4);

%----------coeficientes Ap -  Anb - Scd ---------

Ap = []; %Coeficiente que multiplica a Up, lado izquierdo de la igualda, Salida
Anb_uab =[]; %Producto de los coeficientes Anb y velocidades Ap
S_dc =[]; %Terminos de difusión no ortogonal

%----------Difussion terms---------
%-----Face W -----
%distance between nodes  
%delta_ep_w=sqrt((nc_w(1)-nc_p(1))^2 + (nc_w(2)-nc_p(2))^2);
%non ortogonal diffusion term
s_cd_w= nu*tan_w*0.25*(u_nw+ u_n -u_sw- u_s);
%diffusion term Gradient P_Ai
d_w = nu*len_f(1)/(cos_w*delta_ep_w);

Ap(1)=d_w;
Anb_uab(1)=d_w*u_w;
S_dc(1)=s_cd_w;

%----- Face N -----
delta_ep_n=sqrt((nc_n(1)-nc_p(1))^2 + (nc_n(2)-nc_p(2))^2);
%non ortogonal diffusion term
s_cd_n= nu*tan_n*0.25*(u_ne+ u_e -u_nw- u_w);
%diffusion term Gradient P_Ai
d_n = nu*len_f(2)/(cos_n*delta_ep_n);

Ap(2)=d_n;
Anb_uab(2)=d_n*u_n;
S_dc(1)=s_cd_n;

%----- Face E -----
delta_ep_e=sqrt((nc_e(1)-nc_p(1))^2 + (nc_e(2)-nc_p(2))^2);
%non ortogonal diffusion term
s_cd_e= nu*tan_e*0.25*(u_se+ u_s -u_ne- u_n);
%diffusion term Gradient P_Ai
d_e = nu*len_f(3)/(cos_e*delta_ep_e);

Ap(3)=d_e;
Anb_uab(3)=d_e*u_e;
S_dc(3)=s_cd_e;

%----- Face S -----
delta_ep_s=sqrt((nc_s(1)-nc_p(1))^2 + (nc_s(2)-nc_p(2))^2);
%non ortogonal diffusion term
s_cd_s= nu*tan_s*0.25*(u_sw+ u_w -u_se- u_e);
%diffusion term Gradient P_Ai
d_s = nu*len_f(4)/(cos_s*delta_ep_s);

Ap(4)=d_s;
Anb_uab(4)=d_s*u_s;
S_dc(4)=s_cd_s;

%--------convective terms----------
%mass fluxes for U vel
%face w
fu_w=0.5*(u_w+u_p)*len_f(1);
[Ap,Anb_uab]=esq_interp_upwind(fu_w,u_w,Ap,Anb_uab,fu_w);
%face n
fu_n =0.5*(u_n + u_p)*angles_fns(1,3)*len_f(2)*rho;
%face e
fu_e = -0.5*(u_p + u_e)*len_f(3)*rho;
[Ap,Anb_uab]=esq_interp_upwind(fu_e,u_e,Ap,Anb_uab,fu_e);
%Face S
fu_s =-0.5*(u_p + u_s)*len_f(4)*angles_fns(2,3);

%mass fluxes for v
%face n
fv_n = -0.5*(v_n+v_p)*len_f(2)*angles_fns(1,2)*rho;
[Ap,Anb_uab]=esq_interp_upwind(fv_n,u_n,Ap,Anb_uab,fv_n);
[Ap,Anb_uab]=esq_interp_upwind(fv_n,u_n,Ap,Anb_uab,fu_n);

%face S
fv_s= 0.5*(v_s+v_p)*len_f(4)*angles_fns(2,2)*rho;
[Ap,Anb_uab]=esq_interp_upwind(fv_s,u_s,Ap,Anb_uab,fv_s);
[Ap,Anb_uab]=esq_interp_upwind(fv_s,u_s,Ap,Anb_uab,fu_s);

%Derivate of pressure X
%create vector of diference of transported quantity
delta_P = [p_w;p_n;p_e;p_s] - p_p;
grad_P= wl_op*delta_P;
dpx=grad_P(1);
dpy=grad_P(2);

%--------    Final Equation    ----------
Ap = Ap(~isnan(Ap));
Anb_uab = Anb_uab(~isnan(Anb_uab));
S_dc = S_dc(~isnan(S_dc));

a_P_i = sum(Ap); % Coefficient needed to calculate U_p and pressure correction;
d_Ij = -delt_v/ (a_P_i * dx);
u_star = (1 / a_P_i) * (sum(Anb_uab) + sum(S_dc) - dpx * delt_v); % New U_P, main output
%end





