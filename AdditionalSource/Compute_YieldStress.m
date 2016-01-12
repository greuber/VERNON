function [ NUM ] = Compute_YieldStress(NUM,p,MESH,i,j,PAR )
%% -------------- %% Yield stress computation function %% -------------- %%
% Computes the yield stress by:
% Tau_yield = p*sin(phi)+Cohesion*cos(phi)
%%-----------------------------------------------------------------------%% 

if strcmp(NUM.NUMERICS.Elementtype,'linear') == 1   
else
    p = mean(p);
end

sin_phi = sin(deg2rad(MESH.CompVar.Phi(NUM.Number.number_quad(j,i))));
cos_phi = cos(deg2rad(MESH.CompVar.Phi(NUM.Number.number_quad(j,i))));

Tau_yield = p .* sin_phi + MESH.CompVar.Cohesion(NUM.Number.number_quad(j,i)) .* cos_phi;    % should actually be the mean pressure in my opinion
Tau_yield = max(0.01 * MESH.CompVar.Cohesion(NUM.Number.number_quad(j,i)), Tau_yield);

% store it
Tau_yield(Tau_yield>PAR.Max_yield) = PAR.Max_yield;
NUM.Plasticity.Tau_yield_2d(j,i)   = Tau_yield;

end

