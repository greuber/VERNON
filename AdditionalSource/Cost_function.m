function [ MESH,NUM,PAR,fcost ] = Cost_function( MESH,NUM,PAR )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fcost = (1/2) * ((NUM.r(2:2:NUM.no_nodes*(NUM.ndof-1)) - NUM.r_ini(2:2:NUM.no_nodes*(NUM.ndof-1)))' * (NUM.r(2:2:NUM.no_nodes*(NUM.ndof-1)) - NUM.r_ini(2:2:NUM.no_nodes*(NUM.ndof-1))));

end

