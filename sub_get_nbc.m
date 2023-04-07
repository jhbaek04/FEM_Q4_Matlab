function [ nbc_connectivity ]  =  sub_get_nbc ( x_node , nbc_position )
% The nodes on the natural boundaries are identified, and boundary elements
% are constructed based on the detected nodes.

num_nbc = length ( nbc_position ) ;
nbc_connectivity = cell ( 1 , num_nbc );
for i = 1 : num_nbc
    idx = nbc_position{i} ( x_node(:,1) , x_node(:,2) ) ;
    
    x_nbc = x_node(idx,:);
    if max(x_nbc(:,1)) > min(x_nbc(:,1)) + 1e-8
        [ ~ , ia ]  =  sort ( x_nbc(:,1) );
    else
        [ ~ , ia ]  =  sort ( x_nbc(:,2) );
    end
    idx  =  idx ( ia );
    
    nbc_connectivity{i} = [idx(1:length(idx)-1) , idx(2:length(idx))]';   % The global node indices of the 1st and 2nd nodes of the boundary element.
    % Each row in nbc_connectivity{i} is associated with each boundary
    % element. The 1st column and 2nd column are associated with the 1st
    % node and the 2nd node of the boundary element, respectively.
end

end