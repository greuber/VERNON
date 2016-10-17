function [NUM ] = Compute_ShapeFunctionsVector( COORD,NUM )
%% ---------- %% Shape function preprocessing function %% ---------- %%
% Preassembles the shape functions in a big vector
%%-------------------------------------------------------------------%% 

for j = 1:NUM.NUMERICS.no_nodes_ele
    % Velocity shape functions
    [ NI,dNI ] = get_shape_functions( NUM.NUMERICS.no_nodes_ele,j,COORD );
    NUM.Shape.NI_vector{j} = NI;
    NUM.Shape.dNI_vector{j} = dNI;
    
    % Pressure shape functions
    [ pNI,pdNI ] = get_shape_functions( NUM.NUMERICS.no_nodes_ele_linear,j,COORD );
    NUM.Shape.pNI_vector{j} = pNI;
    NUM.Shape.pdNI_vector{j} = pdNI;
end


% Bliblablubb das hier ist der master 