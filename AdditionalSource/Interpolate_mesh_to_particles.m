function [PARTICLES,MESH,PAR,NUM] = Interpolate_mesh_to_particles(PARTICLES,MESH,interpol_create,PAR,interpol,NUM)


switch interpol
    case 'MtoP'
        rho_interpol      = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.rho,PARTICLES.x,PARTICLES.z,'nearest');
        mu_ref_interpol   = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.mu_ref,PARTICLES.x,PARTICLES.z,'nearest');
        str_ref_interpol  = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.str_ref,PARTICLES.x,PARTICLES.z,'nearest');
        powerlaw_interpol = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.powerlaw,PARTICLES.x,PARTICLES.z,'nearest');
        shmo_interpol     = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.G,PARTICLES.x,PARTICLES.z,'nearest');
        if NUM.Plasticity.Plasticity
            Phi_interpol = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.Phi,PARTICLES.x,PARTICLES.z,'nearest');
            Cohesion_interpol   = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),MESH.CompVar.Cohesion,PARTICLES.x,PARTICLES.z,'nearest');
        end
        
        
        if isfield(NUM,'Result')
            % interpolate velocities
            PARTICLES.Vx = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),NUM.Result.Ux_vec,PARTICLES.x,PARTICLES.z,'nearest');
            PARTICLES.Vz = griddata(MESH.GCOORD(1,:),MESH.GCOORD(2,:),NUM.Result.Uy_vec,PARTICLES.x,PARTICLES.z,'nearest');
        end
        
        if interpol_create
            PARTICLES     = setfield(PARTICLES,'Var',[]);
            PARTICLES.Var = setfield(PARTICLES.Var,'rho',rho_interpol);
            PARTICLES.Var = setfield(PARTICLES.Var,'str_ref',str_ref_interpol);
            PARTICLES.Var = setfield(PARTICLES.Var,'mu_ref',mu_ref_interpol);
            PARTICLES.Var = setfield(PARTICLES.Var,'powerlaw',powerlaw_interpol);
            PARTICLES.Var = setfield(PARTICLES.Var,'G',shmo_interpol);
            if NUM.Plasticity.Plasticity
                PARTICLES.Var = setfield(PARTICLES.Var,'Phi',Phi_interpol);
                PARTICLES.Var = setfield(PARTICLES.Var,'Cohesion',Cohesion_interpol);
            end
        else
            PARTICLES.Var.rho      = rho_interpol;
            PARTICLES.Var.str_ref  = str_ref_interpol;
            PARTICLES.Var.mu_ref   = mu_ref_interpol;
            PARTICLES.Var.powerlaw = powerlaw_interpol;
            PARTICLES.Var.G        = shmo_interpol;
            if NUM.Plasticity.Plasticity
                PARTICLES.Var.Phi      = Phi_interpol;
                PARTICLES.Var.Cohesion = Cohesion_interpol;
            end
            
            % check for values that need to be deleted
            ind_x = find(isnan(PARTICLES.Vx) == 1);
            ind_z = find(isnan(PARTICLES.Vz) == 1);
            ind   = unique([ind_x,ind_z]);
            
            PARTICLES.Vx(ind)     = [];
            PARTICLES.Vz(ind)     = [];
            PARTICLES.x(ind)      = [];
            PARTICLES.z(ind)      = [];
            PARTICLES.phases(ind) = [];
            
            names = fieldnames(PARTICLES.Var);
            % remove particles outside of the domain
            for field  = 1:length(names);
                data            =   getfield(PARTICLES.Var,names{field});
                data(ind)       = [];
                PARTICLES.Var   = setfield(PARTICLES.Var,names{field},data);
            end
        end
        
    case 'PtoM'
        
        rho_interpol      = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.rho,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
        mu_ref_interpol   = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.mu_ref,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
        str_ref_interpol  = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.str_ref,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
        powerlaw_interpol = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.powerlaw,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
        shmo_interpol     = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.G,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
        if NUM.Plasticity.Plasticity
            Phi_interpol      = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.Phi,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
            Cohesion_interpol = griddata(PARTICLES.x,PARTICLES.z,PARTICLES.Var.Cohesion,MESH.GCOORD(1,:),MESH.GCOORD(2,:),'nearest');
        end
        
        
        MESH.CompVar.rho      = rho_interpol;
        MESH.CompVar.str_ref  = str_ref_interpol;
        MESH.CompVar.mu_ref   = mu_ref_interpol;
        MESH.CompVar.powerlaw = powerlaw_interpol;
        MESH.CompVar.G        = shmo_interpol;
        if NUM.Plasticity.Plasticity
            MESH.CompVar.Phi        = Phi_interpol;
            MESH.CompVar.Cohesion   = Cohesion_interpol;
        end
        
        
end
