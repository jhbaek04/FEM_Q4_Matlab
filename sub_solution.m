function [ d ] = sub_solution ( K , f , BC )

%% nodal coefficient vector
d  =  zeros ( size ( f ) ) ;

%% EBC treatment
ebcx_node_index   =  BC.ebcx_node;
ebcy_node_index   =  BC.ebcy_node;

% strong imposition of displacement BC
d ( 2*ebcx_node_index-1 )  =  BC.ebcx_value ;
d ( 2*ebcy_node_index   )  =  BC.ebcy_value ;

% corresponding K indices
idx_ebc   =  sort ( [2*ebcx_node_index-1 ; 2*ebcy_node_index] );
% entire K indices
idx_all   =  [1:size(f)]';
% free indices
idx_free  =  setdiff ( idx_all , idx_ebc ) ;

% reduced K
K_reduced  =  K ( idx_free , idx_free ) ;
% reduced f
f_reduced  =  f ( idx_free )  -  K ( idx_free , idx_ebc ) * d (idx_ebc) ;

%% Solve for d(idx_free)
d ( idx_free )  =  K_reduced \ f_reduced ;


end