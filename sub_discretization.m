function [Mesh] = sub_discretization(Model)

Mesh = struct;

Mesh.x_node = node_distribution ( Model ) ;

Mesh.connectivity = node_connectivity ( Model , Mesh.x_node ) ;

% plot the mesh to check if meshing was correctly done.
plot_mesh ( Mesh.x_node , Mesh.connectivity ); title('Undeformed Mesh')

end


%% nodal distribution
function [x_node] = node_distribution ( Model )
% Nodes are uniformlyl populated.
% The nodal coordinates are stored in x_node.
% x_node : num_nodes - by - 2

x_min = Model.domain.minmax(1,1);
x_max = Model.domain.minmax(2,1);
y_min = Model.domain.minmax(1,2);
y_max = Model.domain.minmax(2,2);

num_elem_x  =  Model.mesh.interval(1);
num_elem_y  =  Model.mesh.interval(2);

x_line     = linspace ( x_min , x_max , num_elem_x+1 )';
y_line     = linspace ( y_min , y_max , num_elem_y+1 )';
one_vector = ones ( size(x_line) );

x_node = [];
for j = 1 : num_elem_y+1
    x_node = [ x_node ; [x_line, y_line(j)*one_vector] ];
end

end

%% element connectivity
function [connectivity] = node_connectivity ( Model, x_node )
% Finite elements are generated, i.e., what global nodes are connected
% to each element is specificed.
% The following numbering rule is used in this code.
% 4 --------- 3
%   |       |
%   |       |
%   |       |
% 1 --------- 2
% The square shown above represents an element, and the numbers denote the
% local node indices. The local nodes are numbered counter-clockwise.
% Then, we will find the global indices associated with each local indices
% and store the information.
% As we are considering a uniform discretization, a simple approach is
% used.


num_elem_x  =  Model.mesh.interval(1) ;
num_elem_y  =  Model.mesh.interval(2) ;
num_element =  num_elem_x * num_elem_y;  % the number of finite elements in the domain

NEN = 4;  % num of nodes per element
connectivity = zeros ( NEN , num_element );
% The i-th row of j-th column in connectivity carries the global node index
% of the local node i of the element j.


% We first identify how big the global node index associated with each 
% local node index is, compared to the global node index of the local node 1.
local_to_global_block = [ 0 ; 1 ; 1+(num_elem_x+1) ; 0+(num_elem_x+1) ] ;
% This only works for the node numbering rule used in the function "node_distribution" defined above.


x_max = Model.domain.minmax(2,1);
y_max = Model.domain.minmax(2,2);

element_index = 0;
for node_index = 1 : length(x_node)  % loop over nodes
    x = x_node(node_index,:);
    if x(1) < x_max-1e-8  &&  x(2) < y_max-1e-8  % skip the nodes on the right and upper boundaries.
        element_index = element_index + 1;
        % The (node-index)-th node is the lower-left-corner node of the
        % currunt element.
        connectivity( : , element_index )  =  node_index + local_to_global_block;
    end
end


end






