function [ N , B ]  =  sub_get_N_and_B ( N0 , dN_dxi , x_node_local )


%% N matrix
N = zeros ( 2 , 8 ) ;
N ( 1 , 1:2:7 )  =  N0;
N ( 2 , 2:2:8 )  =  N0;


%% B matrix
% push forward
F = dN_dxi * x_node_local ; % based on iso-parametric mapping
dN_dx = F \ dN_dxi ;  % chain rule

% B matrix
B = zeros ( 3 , 8 ) ;
B ( 1 , 1:2:7 )  =  dN_dx ( 1 , : ) ;  % dN/dx
B ( 2 , 2:2:8 )  =  dN_dx ( 2 , : ) ;  % dN/dy
B ( 3 , 1:2:7 )  =  dN_dx ( 2 , : ) ;  % dN/dy
B ( 3 , 2:2:8 )  =  dN_dx ( 1 , : ) ;  % dN/dx




end