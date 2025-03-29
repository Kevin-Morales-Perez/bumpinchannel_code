function [xc,yc,cell_face,trig_cell,center_face,cell_volume] = cell_collocated_node(x1,y1)

%this function calculates all geometric parameters of each cell
%funtion to calculate the centroid of each cell
%the cell is divided in 2 triangles and one rectangle, then is
%used the formula xc=sum(xi*ai)/a_total

dx=x1(4)-x1(1);%step in X
dy1= y1(2)-y1(1);%lenght of left face
dy2=y1(3)-y1(4);%lenght of right face
if y1(3)>y1(2)
    a=y1(4)-y1(1);%height of triangle 1
    b=y1(3)-y1(2);%height of triangle 2
    c= y1(2)-y1(4);% height of rectangle 1 

    xt1=(1/3)*dx; % x position of centroid for triangle 1
    xt2=(2/3)*dx; % x position of centroid for triangle 2
    y_in=y1(1);%y origin
    
else
    a=y1(1)-y1(4);%height of triangle 1
    b=y1(2)-y1(3);%height of triangle 2
    c=y1(3)-y1(1);% height of rectangle 1 

    xt1=(2/3)*dx; % x position of centroid for triangle 1
    xt2=(1/3)*dx; % x position of centroid for triangle 2
    y_in=y1(4);%y origin
end

xt3=0.5*dx;% x position of centroid for rectangle 2 

yt1=(2/3)*a;
yt2=a + c + (1/3)*b;
yt3= a + 0.5*c;

at1=0.5*a*dx;%area of triangle 1;
at2=0.5*b*dx;%area of triangle 2;
at3=dx*c;%area of rectangle 3;

a_total=0.5*(dy1+dy2)*dx;%total area of the cell
xc=x1(1)+(xt1*at1 + xt2*at2 + xt3*at3)/a_total;%x centroid of the cell
yc=y_in+ (at1*yt1 + at2*yt2 + at3*yt3)/a_total;%y centroid of the cell

cell_face=[0 0 0 0];%vector which contain lenght of the faces
cell_face(1)=dy1;
cell_face(2)=sqrt((y1(3)-y1(2))^2 + dx^2);
cell_face(3)=dy2;
cell_face(4)=sqrt((y1(4)-y1(1))^2 + dx^2);


%cell volume
cell_volume =0.5*(dy1+dy2)*dx;
%trigonometric funtions of angles from face 2) and 4)
trig_cell=zeros(2, 3);

%face 2
tan1=(y1(3)-y1(2))/dx;
phi1=atan(tan1);
cos1=cos(phi1);
%cos1=1/(1 + tan1^2);
sin1=sin(phi1);
%sin1=tan1/(1 + tan1^2);


trig_cell(1,1)=tan1;
trig_cell(1,2)=cos1;
trig_cell(1,3)=sin1;

%face 4
tan2=(y1(4)-y1(1))/dx;
phi2=atan(tan2);
%cos2=1/(1 + tan2^2);
cos2=cos(phi2);
%sin2=tan2/(1 + tan2^2);
sin2=(phi2);

trig_cell(2,1)=tan2;
trig_cell(2,2)=cos2;
trig_cell(2,3)=sin2;

%cetres of faces
center_face=zeros(4,2);
%face 1
center_face(1,1)=x1(1);
center_face(1,2)=y1(1) + 0.5*dy1;

%face 2
center_face(2,1)=x1(1)+ 0.5*cell_face(2)*cos1;
center_face(2,2)=y1(2) + 0.5*cell_face(2)*sin1;

%face 3
center_face(3,1)=x1(4);
center_face(3,2)=y1(4) + 0.5*dy2;

%face 4
center_face(4,1)=x1(1) + 0.5*cell_face(4)*cos2;
center_face(4,2)=y1(1) + 0.5*cell_face(4)*sin2;
