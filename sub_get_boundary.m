function [BC] = sub_get_boundary ( Model , Mesh )

bc = Model.bc;
x_node = Mesh.x_node;

BC = struct;

%% Essential BC
% Identify the essential boundary nodes and associated displcement values
[ BC.ebcx_node , BC.ebcx_value ] = sub_get_ebc ( x_node , bc.ebcx_position , bc.ebcx_value ) ;
[ BC.ebcy_node , BC.ebcy_value ] = sub_get_ebc ( x_node , bc.ebcy_position , bc.ebcy_value ) ;

%% Natural BC
% Identify the natural boundary nodes and construct the boundary elements
% connectivity.
[ BC.nbc_connectivity ]  =  sub_get_nbc ( x_node , bc.nbc_position );


end