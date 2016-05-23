function [ drdx_temp ] = Adjoint_Res( NUM,PAR,MESH,par,CHAR)
% Computes the derivative of the residual versus design variable with a FD
% approximation

h_max = 1e-10;
h = max([(1e-6*abs(NUM.Adjoint.m(par,1))), h_max]);
index = NUM.Adjoint.index{par};
field = NUM.Adjoint.fields{par};    % the perturbed field in the MESH.CompVar structure

% if strcmp(field,'mu_ref') == 1
%     h = 1e-1*abs(NUM.Adjoint.m(par,1));
% end
if strcmp(field,'Hi') == 1
    h = 0.01;
    MESH_ini = MESH;
end

if isfield(MESH.CompVar,field)
    input_ini = MESH.CompVar.(field);
    input = MESH.CompVar.(field);
elseif isfield(PAR,field)
    input_ini = PAR.(field);
    input = PAR.(field);
else
    display('Input design variable not defined in the residual function')
end

for id1 = 1:2
    NUM.Solve.f_res_brutal = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);
    if id1 == 1
        if isfield(MESH.CompVar,field)
            MESH.CompVar.(field)(unique(index)) = input(unique(index))+h;
        elseif strcmp(field,'rad') == 1
            
            MESH.CompVar.rho         = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(1,1);
            MESH.CompVar.mu_ref      = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(4,1);
            MESH.CompVar.str_ref     = ones(1,NUM.NUMERICS.no_nodes) * (-NUM.ebg)/(1/CHAR.Time);
            MESH.CompVar.powerlaw    = ones(1,NUM.NUMERICS.no_nodes) * PAR.n1;
            MESH.CompVar.Phase       = ones(1,NUM.NUMERICS.no_nodes);
            MESH.CompVar.G           = ones(size(MESH.GCOORD(2,:)))*PAR.G/CHAR.Stress;
            rad = input + NUM.NUMERICS.dx + h;    % in this case we add the minimum dx
            ind = find((MESH.GCOORD(1,:) - (PAR.W/2)).^2 + (MESH.GCOORD(2,:) - (PAR.H/2)).^2 < rad^2);
            MESH.CompVar.rho(ind)         = NUM.Adjoint.m(3,1);
            MESH.CompVar.mu_ref(ind)      = PAR.mu_2/CHAR.Viscosity;
            MESH.CompVar.str_ref(ind)     = (-NUM.ebg)/(1/CHAR.Time);
            MESH.CompVar.powerlaw(ind)    = PAR.n1;
            MESH.CompVar.Phase(ind)       = 2;
            MESH.CompVar.G(ind)           = PAR.G/CHAR.Stress;
            
        elseif strcmp(field,'Hi') == 1
            
            ind_interface = round(PAR.H_interface/NUM.NUMERICS.dz)+1;
            
            PAR.H_interface = input + h;
            
            [n,m] = size(NUM.Number.number_2d);
            a = linspace(0,1,ind_interface);
            ab = linspace(0,1,floor(n-ind_interface));
            b = linspace(ab(end-1),0,floor(n-ind_interface));
            factor = [a b];
            ki = 1;
            kj = m;
            
            for ku = 1:n;
                MESH.GCOORD(2,ki:kj)   = MESH.GCOORD(2,ki:kj) + factor(ku)*h;
                ki = kj+1;
                kj = ku*m;
            end
            
            [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
            
%         elseif strcmp(field,'mu_ref') == 1
%             MESH.CompVar.(field)(unique(index)) = input(unique(index))+h;
%             [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
            
        elseif isfield(PAR,NUM.Adjoint.fields{par}) == 1
            PAR.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = input(unique(index))+h;
        else
            display('Design variable seems to not be defined for the residual function update')
        end
    else
        if isfield(MESH.CompVar,field)
            MESH.CompVar.(field)(unique(index)) = input_ini(unique(index));
        elseif strcmp(field,'rad') == 1
            
            MESH.CompVar.rho         = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(3,1);
            MESH.CompVar.mu_ref      = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(4,1);
            MESH.CompVar.str_ref     = ones(1,NUM.NUMERICS.no_nodes) * (-NUM.ebg)/(1/CHAR.Time);
            MESH.CompVar.powerlaw    = ones(1,NUM.NUMERICS.no_nodes) * PAR.n1;
            MESH.CompVar.Phase       = ones(1,NUM.NUMERICS.no_nodes);
            MESH.CompVar.G           = ones(size(MESH.GCOORD(2,:)))*PAR.G/CHAR.Stress;
            rad = input_ini;
            ind = find((MESH.GCOORD(1,:) - (PAR.W/2)).^2 + (MESH.GCOORD(2,:) - (PAR.H/2)).^2 < rad^2);
            MESH.CompVar.rho(ind)         = NUM.Adjoint.m(1,1);
            MESH.CompVar.mu_ref(ind)      = PAR.mu_2/CHAR.Viscosity ;
            MESH.CompVar.str_ref(ind)     = (-NUM.ebg)/(1/CHAR.Time);
            MESH.CompVar.powerlaw(ind)    = PAR.n1;
            MESH.CompVar.Phase(ind)       = 2;
            MESH.CompVar.G(ind)           = PAR.G/CHAR.Stress;
            
        elseif strcmp(field,'Hi') == 1
            
            PAR.H_interface = input_ini;
            
            MESH = MESH_ini;
            
            [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
            
%         elseif strcmp(field,'mu_ref') == 1
%             MESH.CompVar.(field)(unique(index)) = input_ini(unique(index));
%             [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
            
        elseif isfield(PAR,NUM.Adjoint.fields{par}) == 1
            PAR.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = input_ini(unique(index));
        else
            display('Design variable seems to not be defined for the residual function update')
        end
    end
    
    for i = 1:NUM.NUMERICS.no_elems_global     % compute residual
        f_quad_val = zeros(NUM.NUMERICS.no_nodes_ele*2,1);
        f_line_val = zeros(NUM.NUMERICS.no_nodes_ele_linear,1);
        u = NUM.Solve.r(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1);
        p = NUM.Solve.r(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1);
        [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
        NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1) + f_quad_val;
        NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1) + f_line_val;
    end
    
    % Add boundaries to the system
    for i = 1:1:length(NUM.Boundary.bcdof)
        NUM.Solve.f_res_brutal(NUM.Boundary.bcdof(i)) = 0;
    end
    
    if id1 == 1
        res_plus = NUM.Solve.f_res_brutal;
    else
        res_normal = NUM.Solve.f_res_brutal;
    end
end
drdx_temp = (res_plus - res_normal)./(h);

end

