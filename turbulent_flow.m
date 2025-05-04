%understanding Reynols  navier stokes averaging
clc
close all 
fprintf("####################################################\n")

%velocity with fluctuations 
n_points=101;
t_turb=10;

u = 0 + 0.1*rand(1,n_points); 
v = 4.5  + 1*rand(1,n_points); 
w = 0 + 0.1*rand(1,n_points); 

t=linspace(0,t_turb,n_points);
zero_vector=zeros(1,n_points);%for reference

figure(1)
plot(t,u)
title("velocities")
hold on
plot(t,v)
hold on
plot(t,w)
hold on
plot(t,zero_vector)
axis equal
hold off

%mean values
u_mean=(1/t_turb)*trapz(t,u);
fprintf("U mean=")
disp(u_mean)
v_mean=(1/t_turb)*trapz(t,v);
fprintf("V mean=")
disp(v_mean)
w_mean=(1/t_turb)*trapz(t,w);
fprintf("W mean=")
disp(w_mean)

%fluctuating values
u_f=u-u_mean;
v_f=v-v_mean;
w_f=w-w_mean;

figure(2)

plot(t,u_f)
title("Fluctuating values")
hold on
plot(t,v_f)
hold on
plot(t,w_f)
hold on
plot(t,zero_vector)
axis equal
hold off

%averages of fluctuations
u_t_mean=(1/t_turb)*trapz(t,u_f);
fprintf("U fluctuating mean=")
disp(u_t_mean)
v_t_mean=(1/t_turb)*trapz(t,v_f);
fprintf("V fluctuating mean=")
disp(v_t_mean)
w_t_mean=(1/t_turb)*trapz(t,w_f);
fprintf("W fluctuating mean=")
disp(w_t_mean)

%variances
u_f2=u_f.^2;
v_f2=v_f.^2;
w_f2=w_f.^2;

u_f2_mean=(1/t_turb)*trapz(t,u_f2);
fprintf("u'2 variance=")
disp(u_f2_mean)
v_f2_mean=(1/t_turb)*trapz(t,v_f2);
fprintf("v'2 variance=")
disp(v_f2_mean)
w_f2_mean=(1/t_turb)*trapz(t,w_f2);
fprintf("w'2 variance=")
disp(w_f2_mean)

%r.m.s
u_f2_rms=sqrt(u_f2_mean);
fprintf("u'2 rms=")
disp(u_f2_rms)
v_f2_rms=sqrt(v_f2_mean);
fprintf("v'2 rms=")
disp(v_f2_rms)
w_f2_rms=sqrt(w_f2_mean);
fprintf("w'2 rms=")
disp(w_f2_rms)

%Kinetic turbulent energy
k_turb=0.5*(u_f2_mean+v_f2_mean+w_f2_mean);
fprintf("Turbulent kinetik energy =")
disp(k_turb)

%second moments
%u'v'
u_v_f=u_f.*v_f;
u_v_mean=(1/t_turb)*trapz(t,u_v_f);
fprintf("u'v'- =")
disp(u_v_mean)

%u'w'
u_w_f=u_f.*w_f;
u_w_mean=(1/t_turb)*trapz(t,u_w_f);
fprintf("u'w'- =")
disp(u_w_mean)

%v,w
v_w_f=v_f.*w_f;
v_w_mean=(1/t_turb)*trapz(t,v_w_f);
fprintf("v'w'- =")
disp(v_w_mean)


%higher order moments

%Skewness
u_f3=u_f.^3;
v_f3=v_f.^3;
w_f3=w_f.^3;

u_f3_mean=(1/t_turb)*trapz(t,u_f3);
fprintf("u'3 skewness=")
disp(u_f3_mean)
v_f3_mean=(1/t_turb)*trapz(t,v_f3);
fprintf("v'3 skewness=")
disp(v_f3_mean)
w_f3_mean=(1/t_turb)*trapz(t,w_f3);
fprintf("w'3 skewness=")
disp(w_f3_mean)



%Kurtosis (Peakedness)
u_f4=u_f.^4;
v_f4=v_f.^4;
w_f4=w_f.^4;

u_f4_mean=(1/t_turb)*trapz(t,u_f4);
fprintf("u'4 kurtosis=")
disp(u_f3_mean)
v_f4_mean=(1/t_turb)*trapz(t,v_f4);
fprintf("v'4 kurtosis=")
disp(v_f3_mean)
w_f4_mean=(1/t_turb)*trapz(t,w_f4);
fprintf("w'4 kurtosis=")
disp(w_f4_mean)


%Turbulent Stress tensor

turb_stress_tensor=[u_f2_mean,u_v_mean,u_w_mean;u_v_mean,v_f2_mean,v_w_mean;u_w_mean,v_w_mean,w_f2_mean];
fprintf("Reynolds stress tensor = \n")
disp(turb_stress_tensor)


%Representing as a ellipsoid
[V,D]=eig(turb_stress_tensor);


% Get semi-axis lengths (sqrt of eigenvalues for energy ellipsoid)
a = sqrt(D(1,1)); % x-axis
b = sqrt(D(2,2)); % y-axis
c = sqrt(D(3,3)); % z-axis

% Generate a unit sphere
[x, y, z] = ellipsoid(0, 0, 0, 1, 1, 1, 50);

% Scale and rotate the sphere into an ellipsoid
X = V(1,1)*a*x + V(1,2)*b*y + V(1,3)*c*z;
Y = V(2,1)*a*x + V(2,2)*b*y + V(2,3)*c*z;
Z = V(3,1)*a*x + V(3,2)*b*y + V(3,3)*c*z;

% Plot the ellipsoid
figure;
surf(X, Y, Z, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
axis equal;
xlabel('u'); ylabel('v'); zlabel('w');
title('Reynolds Stress Tensor Ellipsoid');
colormap('jet');
colorbar;
grid on;
view(3);