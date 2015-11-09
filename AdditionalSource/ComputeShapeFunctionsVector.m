function [NUM ] = ComputeShapeFunctionsVector( COORD,NUM )
% Computes the shape functions as big vector

for j = 1:NUM.NUMERICS.no_nodes_ele
    % Velocity shape functions
    [ NI,dNI ] = Shape_Functions( NUM.NUMERICS.no_nodes_ele,j,COORD );
    NUM.Shape.NI_vector{j} = NI;
    NUM.Shape.dNI_vector{j} = dNI;
    
    % Pressure shape functions
    [ pNI,pdNI ] = Shape_Functions( NUM.NUMERICS.no_nodes_ele_linear,j,COORD );
    NUM.Shape.pNI_vector{j} = pNI;
    NUM.Shape.pdNI_vector{j} = pdNI;
end


