function [ MESH ] = Compute_StrainWeakening( NUM,i,j,MESH )
% adds strain weakening on the full strainrate 

C_end   = MESH.CompVar.C_end(NUM.Number.number_quad(j,i));
phi_end = MESH.CompVar.phi_end(NUM.Number.number_quad(j,i));

C_ini   = MESH.CompVar.C_ini(NUM.Number.number_quad(j,i));
phi_ini = MESH.CompVar.phi_ini(NUM.Number.number_quad(j,i));

Strain_start    = NUM.Plasticity.Weaken.Strain_start;
Strain_end      = NUM.Plasticity.Weaken.Strain_end;

if NUM.Strain.APS(j,i) > NUM.Plasticity.Weaken.Strain_end
    C_fac       = C_end;
    phi_fac     = phi_end;
    
elseif NUM.Strain.APS(j,i) < Strain_end && NUM.Strain.APS(j,i) > Strain_start
    fac         = (NUM.Strain.APS(j,i) - Strain_start)/(Strain_end - Strain_start);
    C_fac       = C_ini + fac   * (C_end - C_ini);
    phi_fac     = phi_ini + fac * (phi_end - phi_ini);
    
else
    C_fac       = C_ini;
    phi_fac     = phi_ini;
    
end


MESH.CompVar.Cohesion(NUM.Number.number_quad(j,i)) = C_fac;
MESH.CompVar.Phi(NUM.Number.number_quad(j,i))      = phi_fac;