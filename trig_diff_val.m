function [trig_diff_terms] =trig_diff_val(node_cord,angles_fns)
%funtion to obtain trigonometric values for diffusion terms in momentum
%equations
%the output has this format [tan_w cos_w;tan_w cos_w;tan_w cos_w;tan_w
%cos_w] 4 rows 2 columns

%angle phi tangents ( used for difussion eq)
%example of node imputs 

%trigonometric funtions for faces 1 ,N) and 2 ,S)
%angles_fns= [0 0 0;0 0 0];
%nc_w=node_cord(1,:);
%nc_n=node_cord(2,:);
%nc_e=node_cord(3,:);
%nc_s=node_cord(4,:);
%nc_p=node_cord(5,:);

%tan_w=(nc_w(2)-nc_p(2))/(nc_w(1)-nc_p(1));
tan_w=(node_cord(1,2)-node_cord(5,2))/(node_cord(1,1)-node_cord(5,1));
phi_w=atan(tan_w);
cos_w=cos(phi_w);

tan_n=tan(atan(angles_fns(1,1)) + atan2(node_cord(2,1)-node_cord(5,1),node_cord(2,2)-node_cord(5,2)));
phi_n=atan(tan_n);
cos_n=cos(phi_n);

tan_e=(node_cord(2,2)-node_cord(5,2))/node_cord(3,1)-node_cord(5,1);
phi_e=atan(tan_e);
cos_e=cos(phi_e);

tan_s =tan(atan(angles_fns(2,1)) +  atan2(node_cord(5,1)-node_cord(4,1),node_cord(5,2)-node_cord(4,2)));
phi_s=atan(tan_s);
cos_s=cos(phi_s);

trig_diff_terms= [tan_w cos_w;tan_n cos_n;tan_e cos_e;tan_s,cos_s];


