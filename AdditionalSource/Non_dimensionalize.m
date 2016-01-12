function [ NUM,CHAR,PAR,MESH,varargout] = Non_dimensionalize( NUM,CHAR,PAR,MESH ,varargin)
%% -------------- %% Nondimensionalization function %% -------------- %%
% Nondimensionalization for every parameter
%%--------------------------------------------------------------------%% 

if size(varargin)>=1
    PARTICLES = varargin{1};
end

if ~isfield(CHAR,'Length')
    CHAR.Length         =       100e3;
end
if ~isfield(CHAR,'Viscosity')
    CHAR.Viscosity      =       1e20;
end
if ~isfield(CHAR,'Time')
    CHAR.Time           =       1/1e-15;
end

% Derived characteristic values
CHAR.Velocity           =       CHAR.Length/CHAR.Time;
CHAR.Stress             =       CHAR.Viscosity/CHAR.Time;
CHAR.kg                 =       CHAR.Stress*CHAR.Length*CHAR.Time^2;
CHAR.rho                =       CHAR.kg/CHAR.Length^3;
CHAR.Gravity            =       CHAR.Length/CHAR.Time^2;

if isfield(NUM ,'ebg')
    NUM.ebg = NUM.ebg/(1/CHAR.Time);
end

if isfield(NUM ,'vbg')
    NUM.vbg = NUM.vbg/CHAR.Velocity;
end

if isfield(MESH ,'GCOORD')
    MESH.GCOORD = MESH.GCOORD./CHAR.Length;
end

if isfield(MESH.CompVar, 'mu_ref')
    MESH.CompVar.mu_ref = MESH.CompVar.mu_ref./CHAR.Viscosity;         % prefactor viscosity
end

if isfield(MESH.CompVar,'str_ref')
    MESH.CompVar.str_ref  = MESH.CompVar.str_ref./(1/CHAR.Time);           % characteristic strainrate
end

if isfield(MESH.CompVar, 'rho')
    MESH.CompVar.rho      = MESH.CompVar.rho./CHAR.rho;
end

if isfield(MESH.CompVar, 'Cohesion')
    MESH.CompVar.Cohesion      = MESH.CompVar.Cohesion./CHAR.Stress;
end

if isfield(MESH.CompVar, 'C_ini')
    MESH.CompVar.C_ini      = MESH.CompVar.C_ini./CHAR.Stress;
end

if isfield(MESH.CompVar, 'C_end')
    MESH.CompVar.C_end      = MESH.CompVar.C_end./CHAR.Stress;
end

if isfield(MESH.CompVar, 'G')
    MESH.CompVar.G      = MESH.CompVar.G./CHAR.Stress;
end

if isfield(NUM.Viscosity, 'UpperCutoff')
    NUM.Viscosity.UpperCutoff      = NUM.Viscosity.UpperCutoff./CHAR.Viscosity;
end

if isfield(NUM.Viscosity, 'LowerCutoff')
    NUM.Viscosity.LowerCutoff      = NUM.Viscosity.LowerCutoff./CHAR.Viscosity;
end

if isfield(NUM.Plasticity , 'Max_yield')
    NUM.Plasticity.Max_yield       = NUM.Plasticity.Max_yield ./CHAR.Stress;
end


if exist('PARTICLES','var')
    if isfield(PARTICLES.Var, 'mu_ref')
        PARTICLES.Var.mu_ref = PARTICLES.Var.mu_ref./CHAR.Viscosity;         % prefactor viscosity
    end
    
    if isfield(PARTICLES.Var,'str_ref')
        PARTICLES.Var.str_ref  = PARTICLES.Var.str_ref./(1/CHAR.Time);           % characteristic strainrate
    end
    
    if isfield(PARTICLES.Var, 'rho')
        PARTICLES.Var.rho      = PARTICLES.Var.rho./CHAR.rho;
    end
    
    if isfield(PARTICLES.Var, 'Cohesion')
        PARTICLES.Var.Cohesion      = PARTICLES.Var.Cohesion./CHAR.Stress;
    end
    
    if isfield(PARTICLES.Var, 'G')
        PARTICLES.Var.G      = PARTICLES.Var.G./CHAR.Stress;
    end
    
    PARTICLES.x = PARTICLES.x./CHAR.Length;
    PARTICLES.z = PARTICLES.z./CHAR.Length;
    
    varargout{1} = PARTICLES;
end


if isfield(PAR,'g')
    PAR.g = PAR.g/CHAR.Gravity;
end

if isfield(PAR,'Max_yield')
    PAR.Max_yield = PAR.Max_yield/CHAR.Stress;
end

if isfield(PAR, 'dt')
    PAR.dt = PAR.dt/CHAR.Time;
end

if isfield(PAR, 'dt_max')
    PAR.dt_max = PAR.dt_max/CHAR.Time;
end

if isfield(NUM.NUMERICS ,'dz')
    NUM.NUMERICS.dz = NUM.NUMERICS.dz/CHAR.Length;
end

if isfield(NUM.NUMERICS, 'dx')
    NUM.NUMERICS.dx = NUM.NUMERICS.dx/CHAR.Length;
end

if isfield(PAR, 'W')
    PAR.W = PAR.W/CHAR.Length;
end

if isfield(PAR, 'H')
    PAR.H = PAR.H/CHAR.Length;
end

if isfield(PAR, 'A0')
    PAR.A0 = PAR.A0/CHAR.Length;
end



end

