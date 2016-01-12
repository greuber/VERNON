function [ NUM ] = get_numerical_parameters( NUM,PAR )
% Computes the numerical parameters
NUM.NUMERICS.no_elems_global = NUM.NUMERICS.no_elems_x*NUM.NUMERICS.no_elems_z;

if ~isfield(NUM.NUMERICS,'Elementtype')
    NUM.NUMERICS.Elementtype = 'quadratic';
end

switch NUM.NUMERICS.Elementtype
    case 'linear'
        NUM.NUMERICS.no_nodes_ele = 4;
        NUM.NUMERICS.no_nodes_ele_linear = 1;
        NUM.NUMERICS.no_intp = 4;        % 4 or 9
        NUM.NUMERICS.no_intp_linear = 1;
        NUM.NUMERICS.ndof = 3;
        NUM.NUMERICS.Dimension = 2;
        
        NUM.NUMERICS.no_nodes_x = (NUM.NUMERICS.no_elems_x)+1;
        NUM.NUMERICS.no_nodes_z = (NUM.NUMERICS.no_elems_z)+1;
        NUM.NUMERICS.no_nodes = NUM.NUMERICS.no_nodes_x*NUM.NUMERICS.no_nodes_z;
        
        NUM.NUMERICS.no_nodes_x_linear = (NUM.NUMERICS.no_elems_x);
        NUM.NUMERICS.no_nodes_z_linear = (NUM.NUMERICS.no_elems_z);
        
    case 'quadratic'
        NUM.NUMERICS.no_nodes_ele = 9;
        NUM.NUMERICS.no_nodes_ele_linear = 4;
        NUM.NUMERICS.no_intp = 9;        % 4 or 9
        NUM.NUMERICS.no_intp_linear = 4;
        NUM.NUMERICS.ndof = 3;
        NUM.NUMERICS.Dimension = 2;
        
        NUM.NUMERICS.no_nodes_x = (NUM.NUMERICS.no_elems_x*2)+1;
        NUM.NUMERICS.no_nodes_z = (NUM.NUMERICS.no_elems_z*2)+1;
        NUM.NUMERICS.no_nodes = NUM.NUMERICS.no_nodes_x*NUM.NUMERICS.no_nodes_z;
        
        NUM.NUMERICS.no_nodes_x_linear = (NUM.NUMERICS.no_elems_x)+1;
        NUM.NUMERICS.no_nodes_z_linear = (NUM.NUMERICS.no_elems_z)+1;
end

NUM.NUMERICS.dx = (PAR.W/(NUM.NUMERICS.no_nodes_x-1));
NUM.NUMERICS.dz = (PAR.H/(NUM.NUMERICS.no_nodes_z-1));

end

