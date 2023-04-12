function [ N , dN_dxi1 , dN_dxi2 ]  =  sub_shape_local ( xi )
%% input
% xi: npt - by - 2, parametric coordinate of a set of evaluation points.
% npt: num of evaluation points

%% output
% N:    npt - by - nen
% dN1:  npt - by - nen
% dN2:  npt - by - nen
% nen: num of node per element

% N = [ 4-entry shape function at 1st evaluataion point ;
%       4-entry shape function at 2nd evaluataion point ;
%                             :
%       4-entry shape function at npt-th evaluataion point ]

%%

N = .25 * [(1-xi(:,1)).*(1-xi(:,2)), (1+xi(:,1)).*(1-xi(:,2)), (1+xi(:,1)).*(1+xi(:,2)), (1-xi(:,1)).*(1+xi(:,2))];

dN_dxi1 = .25 * [-(1-xi(:,2)), +(1-xi(:,2)), +(1+xi(:,2)), -(1+xi(:,2))];
dN_dxi2 = .25 * [-(1-xi(:,1)), -(1+xi(:,1)), +(1+xi(:,1)), +(1-xi(:,1))];

end
