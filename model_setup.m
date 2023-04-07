function [Model] = model_setup
% Set up the model.

Model = struct;

%% material properties
E       =  100e3;
nu      =  0.3;
lambda  =  E * nu / (1+nu) / (1-2*nu);
mu      =  E / 2 / (1+nu);
Model.material.Youngs_modulus         =  E;
Model.material.Poissons_ratio         =  nu;
Model.material.Lames_first_parameter  =  lambda;
Model.material.shear_modulus          =  mu;

%% geometry of domain
% (x_min, y_max) ----------------- (x_max, y_max)
%                |               |
%                |               |
%                |               |
%                |               |
%                |               |
% (x_min, y_min) ----------------- (x_max, y_min)
x_min = -1; x_max = +1;
y_min = -1; y_max = +1;
Model.domain.minmax = [[x_min; x_max], [y_min; y_max]];

%% mesh info
% for uniform mesh. The interval in x and y directions are specified.
% The actural meshing will be done later in "sub_discretization".
interval_x = 10; interval_y = 10;
Model.mesh.interval = [interval_x, interval_y];

%% Define boundary conditions
% One boundary segment can have different type of BCs (essential or
% natural) in different directions (x or y) in this code.
% The essential BC (EBC) is seperately defined for each direction.
% The natural BC (NBC) should be specified for both directions, but put
% zero for the direction in which an EBC is specified.

%% Essential boundaries - x direction
num_ebcx =  2;  % num of essential boundaries.
Model.bc.ebcx_position    =  cell ( 1 , num_ebcx ) ;
Model.bc.ebcx_value       =  cell ( 1 , num_ebcx ) ;

idx_bc = 0;
% x essential (Dirichlet) boundary 1
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.ebcx_position{idx_bc} =  @(x,y) find( abs(x-x_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.ebcx_value{idx_bc}    =  @(x,y) [ 0.1*x+0.3*y ] .* ones(size(x));  % displacement in x direction

% x essential (Dirichlet) boundary 2
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.ebcx_position{idx_bc} =  @(x,y) find( abs(y-y_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.ebcx_value{idx_bc}    =  @(x,y) [ 0.1*x+0.3*y ] .* ones(size(x));  % displacement in x direction

%% Essential boundaries - y direction
num_ebcy =  2;  % num of essential boundaries.
Model.bc.ebcy_position    =  cell ( 1 , num_ebcy ) ;
Model.bc.ebcy_value       =  cell ( 1 , num_ebcy ) ;

idx_bc = 0;
% y essential (Dirichlet) boundary 1
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.ebcy_position{idx_bc} =  @(x,y) find( abs(x-x_max) < 1e-8 );  % where the condition below will be applied.
Model.bc.ebcy_value{idx_bc}    =  @(x,y) [ 0.2*x+0.4*y ] .* ones(size(x));  % displacement in y direction

% y essential (Dirichlet) boundary 2
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.ebcy_position{idx_bc} =  @(x,y) find( abs(y-y_max) < 1e-8 );  % where the condition below will be applied.
Model.bc.ebcy_value{idx_bc}    =  @(x,y) [ 0.2*x+0.4*y ] .* ones(size(x));  % displacement in y direction

%% Natural boundaries - for both directions
num_nbc    =  4;  % num of natural boundaries.
Model.bc.nbc_position    =  cell ( 1 , num_nbc ) ;
Model.bc.nbc_value       =  cell ( 1 , num_nbc ) ;

idx_bc = 0;
% natural (Neumann) boundary 1
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.nbc_position{idx_bc} = @(x,y) find( abs(x-x_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.nbc_value{idx_bc}    = @(x,y) -[ 0  ,  0.5*mu ] .* ones(size(x));  % traction in x and y directions. Put zero if an essential BC is applied to that direction.

% natural (Neumann) boundary 2
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.nbc_position{idx_bc} = @(x,y) find( abs(y-y_min) < 1e-8 );  % where the condition below will be applied.
Model.bc.nbc_value{idx_bc}    = @(x,y) -[ 0  ,  0.4*(2*mu+lambda)+0.1*lambda ] .* ones(size(x));  % traction in x and y directions. Put zero if an essential BC is applied to that direction.

% natural (Neumann) boundary 3
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.nbc_position{idx_bc} = @(x,y) find( abs(x-x_max) < 1e-8 );  % where the condition below will be applied.
Model.bc.nbc_value{idx_bc}    = @(x,y) [ 0.1*(2*mu+lambda)+0.4*lambda  ,  0 ] .* ones(size(x));  % traction in x and y directions. Put zero if an essential BC is applied to that direction.

% natural (Neumann) boundary 4
% The actual (x,y) values are plugged in after the domain is meshed.
idx_bc = idx_bc + 1;
Model.bc.nbc_position{idx_bc} = @(x,y) find( abs(y-y_max) < 1e-8 );  % where the condition below will be applied.
Model.bc.nbc_value{idx_bc}    = @(x,y) [ 0.5*mu  ,  0 ] .* ones(size(x));  % traction in x and y directions. Put zero if an essential BC is applied to that direction.

% If neighter EBC nor NBC is specified for a boundary segment, that
% boundary segment is treated as a natural boundary with a zero traction.

%% Body force
% The actual (x,y) values are plugged in after the domain is meshed.
Model.body_force   =  @(x,y) [ 0 , 0 ] .* ones(size(x));

%% Exact solution - if there exists
% The actual (x,y) values are plugged in after the domain is meshed.
Model.exact.use         =   1;  % 0 - exact solution not used ,   1 - used
Model.exact.displ       =   @(x,y) [ 0.1*x+0.3*y  ,  0.2*x+0.4*y ] .* ones(size(x));
Model.exact.strain      =   @(x,y) [ 0.1  ,  0.4 ,  0.5 ] .* ones(size(x));
Model.exact.stress      =   @(x,y) [ 0.1*(2*mu+lambda)+0.4*lambda  , ...
                                     0.4*(2*mu+lambda)+0.1*lambda  , ...  
                                     0.5*mu ] .* ones(size(x));

end