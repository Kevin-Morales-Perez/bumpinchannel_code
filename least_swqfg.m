function [geof_lws] = least_swqfg(cent_adj,cent_p)
%Funtion to obtain geometrical parameters needed for weighted least squares
%Gradient
%Centroids from adyacent nodes N, faces enumerates as following [1,2,3,4]
%example input cent_adj=[-1,1;1,3;3,1;1,-2];
%centroid of cell P
%example input cent_p = [1,1];
%Matrix of distance vector from node P to nodes N
dist_pn=cent_adj - ones(4,1)*cent_p;
w_D=zeros(4);
%matrix of weights
for i=1:4
    for j=1:4
        if i==j
            w_D(i,j)=1/norm(dist_pn(i,:));
        end
    end
end

g = (dist_pn'*w_D')*w_D*dist_pn;%g matrix , this are the matrices that multiplies the gradient at the left side of the equation
geof_lws=g^-1*dist_pn'*(w_D'*w_D);%parameter needed to be multiplied per vector of diferences of values to obtain gradient at node P
end