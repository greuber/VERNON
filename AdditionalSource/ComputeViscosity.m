function [ NUM,MESH] = ComputeViscosity(i,NUM,j,MESH,CHAR,PAR,B)
% Computes viscosity

LowerCutoff = NUM.Viscosity.LowerCutoff;
UpperCutoff = NUM.Viscosity.UpperCutoff;

if isfield(NUM.Viscosity,'mu_2d')
else
    NUM.Viscosity.mu_2d         = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Strain.str_invariant_2d = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Strain.E2nd_total       = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Stress.T2nd             = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Stress.T2nd_old         = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Plasticity.Tau_yield_2d = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

if isfield(NUM.Stress,'Txx')
else
    NUM.Stress.Txx          = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Stress.Tzz          = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Stress.Txz          = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

if isfield(NUM.Stress,'Txx_old')
else
    NUM.Stress.Txx_old = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Stress.Tzz_old = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Stress.Txz_old = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

if isfield(NUM.Strain,'Exx')
else
    NUM.Strain.Exx = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Strain.Ezz = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    NUM.Strain.Exz = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

if NUM.Plasticity.Plasticity
    if isfield(NUM.Strain,'APS')
    else
        NUM.Strain.APS = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    end
    
    if isfield(NUM.Strain,'E_pl_2nd')
    else
        NUM.Strain.E_pl_2nd = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    end
    
    if isfield(NUM.Viscosity,'mu_pl_2d')
    else
    NUM.Viscosity.mu_pl_2d = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
    end
end

%% Elastic strainrate

NUM.Stress.tau_local =  [NUM.Stress.Txx_old(j,i), NUM.Stress.Tzz_old(j,i) ,NUM.Stress.Txz_old(j,i)]';

NUM.Strain.str_local   = B*NUM.Solve.r(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i));
NUM.Strain.str_local(3,1) = NUM.Strain.str_local(3,1)/2;
NUM.Strain.str_total = NUM.Strain.str_local;
NUM.Strain.str_local = NUM.Strain.str_local + (NUM.Stress.tau_local/(2*MESH.CompVar.G(NUM.Number.number_quad(j,i))*PAR.dt));

NUM.Strain.str_invariant = sqrt((1/2)*(NUM.Strain.str_local(1,1)^2 + NUM.Strain.str_local(2,1)^2) + NUM.Strain.str_local(3,1)^2);
NUM.Strain.str_total_inv = sqrt((1/2)*(NUM.Strain.str_total(1,1)^2 + NUM.Strain.str_total(2,1)^2) + NUM.Strain.str_total(3,1)^2);

if NUM.Strain.str_invariant == 0
    NUM.Strain.str_invariant = MESH.CompVar.str_ref(NUM.Number.number_quad(j,i));
end

NUM.Strain.Exx(j,i) = NUM.Strain.str_local(1);
NUM.Strain.Ezz(j,i) = NUM.Strain.str_local(2);
NUM.Strain.Exz(j,i) = NUM.Strain.str_local(3);

%% Nonlinear + elastic part of viscosity

%%%%%%%%%%%%%%%%%% NEEDS TO BE TURNED ON %%%%%%%%%%%%
% NUM.Viscosity.mu = MESH.CompVar.mu_ref(NUM.Number.number_quad(j,i))*(NUM.Strain.str_invariant./MESH.CompVar.str_ref(NUM.Number.number_quad(j,i))).^(1./MESH.CompVar.powerlaw(NUM.Number.number_quad(j,i))-1);
% NUM.Viscosity.mu = 1./(1./NUM.Viscosity.mu + 1./(MESH.CompVar.G(NUM.Number.number_quad(j,i))*PAR.dt)); 
%%%%%%%%%%%%%%%%%% NEEDS TO BE TURNED ON %%%%%%%%%%%%

% For the visco-elastic case uncomment this
NUM.Viscosity.mu = MESH.CompVar.G(NUM.Number.number_quad(j,i)) *PAR.dt;          

mu_vis = NUM.Viscosity.mu;

%% Plastic
if NUM.Plasticity.Plasticity
    
    if NUM.Solve.number_pic == 1 && NUM.time == 0   % in the first iteration in the first timestep pressure is set to 1
        p = ones(size(NUM.Solve.r(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1)));
    else
        p = NUM.Solve.r(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1);
    end
    
    if isfield(NUM.Plasticity,'Weaken')
        [ MESH ] = Compute_StrainWeakening( NUM,i,j,MESH );
    end
    
    [ NUM ] = Compute_YieldStress( NUM,p,MESH,i,j,PAR );
    
    mu_pl  = NUM.Plasticity.Tau_yield_2d(j,i)/(2*NUM.Strain.str_invariant);
    NUM.Viscosity.mu = min(mu_vis, mu_pl);
    
    if NUM.Viscosity.mu == mu_pl
        NUM.Plasticity.Plastic(j,i) = 1;
    end
    
    NUM.Viscosity.mu_pl_2d(j,i) = mu_pl;
    
end

%% Stresses and bounds
% Check boundaries
NUM.Viscosity.mu(NUM.Viscosity.mu<LowerCutoff) = LowerCutoff;
NUM.Viscosity.mu(NUM.Viscosity.mu>UpperCutoff) = UpperCutoff;
NUM.Viscosity.mu(isnan(NUM.Viscosity.mu))      = UpperCutoff;

NUM.Stress.Txx(j,i) = 2*NUM.Viscosity.mu * NUM.Strain.str_local(1,1);
NUM.Stress.Tzz(j,i) = 2*NUM.Viscosity.mu * NUM.Strain.str_local(2,1);
NUM.Stress.Txz(j,i) = 2*NUM.Viscosity.mu * NUM.Strain.str_local(3,1);
NUM.Stress.T2nd(j,i) = sqrt((1/2)*(NUM.Stress.Txx(j,i)^2 + NUM.Stress.Tzz(j,i)^2) + NUM.Stress.Txz(j,i)^2);

NUM.Viscosity.mu_2d(j,i)         = NUM.Viscosity.mu;
NUM.Viscosity.mu_vis_2d(j,i)     = mu_vis;
NUM.Strain.str_invariant_2d(j,i) = NUM.Strain.str_invariant;
NUM.Strain.E2nd_total(j,i)       = NUM.Strain.str_total_inv;    

end
