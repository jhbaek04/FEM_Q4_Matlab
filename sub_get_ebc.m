function [ ebc_node , ebc_value ] = sub_get_ebc ( x_node , ebc_position , ebc_value ) 

num_ebc = length ( ebc_position ) ;
idx = [];
val = [];
for i = 1 : num_ebc
    idx0 = ebc_position{i}( x_node(:,1) , x_node(:,2) );
    idx = [ idx ;  idx0 ] ;
    val = [ val ; ebc_value{i}( x_node(idx0,1) , x_node(idx0,2) ) ] ;
end
[~, ia, ~] = unique(idx);  % corner nodes are shared by multiple boundary segments, so unify the duplicated nodes.
ebc_node = idx(ia);   % Global indices of the nodes that are on the essential boundaries.
ebc_value = val(ia, :);   % corresponding displacement values.


end