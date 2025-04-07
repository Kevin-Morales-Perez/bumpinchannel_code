function [gradient_weights_op,norm] =gradient_lsq_weights(adj_nodes,p_node)
%Function that returns least square gradient operator
%also return the norms of the vectors that join the nodes 
%Centroids from adyacent nodes N, faces enumerates as following
%[1,2,3,4][X,Y]
%example input adj_nodes=[-1,1;1,3;3,1;1,-2];
%centroid of cell P
%example input p_node = [1,1];
%Matrix of distance vector from node P to nodes N
dist_pn=adj_nodes - p_node;
%Vector of weights
norm=vecnorm(dist_pn,2,2);
w_v=1./norm;
w_v=diag(w_v);
%g matrix , this are the matrices that multiplies the gradient at the left side of the equation
G_m=dist_pn'*(w_v*dist_pn);
gradient_weights_op=G_m\(dist_pn'*w_v);
end

