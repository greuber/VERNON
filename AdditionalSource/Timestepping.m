function [ NUM,MESH] = Timestepping(NUM,PAR,MESH,CHAR )

MESH_Old = MESH;

switch NUM.Timestep.method
    
    case 'EulerImplicit'
        
        [ NUM,MESH] = Solver(NUM,PAR,MESH,CHAR );
    
    case 'Euler'
        
        [ NUM,MESH] = Solver(NUM,PAR,MESH,CHAR );
        
    case 'RK2'
        Vx_runge    = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-2),2);
        Vz_runge    = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-2),2);
        P_runge     = zeros(NUM.NUMERICS.no_nodes_linear,2);
        PAR.dt_ini  = PAR.dt;
        
        for Runge = 1:1:2
            if Runge == 1
                PAR.dt = PAR.dt_ini/2;
            else
                PAR.dt = PAR.dt_ini;
            end
            [ NUM,MESH] = Solver(NUM,PAR,MESH,CHAR );
            
            Vx_runge(:,Runge) = NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes));
            Vz_runge(:,Runge) = NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes));
            P_runge(:,Runge)  = NUM.Solve.r(NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+1:end));
            
            MESH.GCOORD(1,:) = MESH.Old(1,:) + Vx_runge(:,Runge)'*PAR.dt;
            MESH.GCOORD(2,:) = MESH.Old(2,:) + Vz_runge(:,Runge)'*PAR.dt;
            
            display(['----- RK2 step ',num2str(Runge),' finished -----']) 
            
        end
        
        MESH = MESH_Old;
        
        PAR.dt = PAR.dt_ini;
        NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes))      = Vx_runge(:,2)';   % Vx update
        NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes))      = Vz_runge(:,2)';
        NUM.Solve.r(NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+1:end))  = P_runge(:,2);
        
    case 'RK4'
        Vx_runge    = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-2),4);
        Vz_runge    = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-2),4);
        P_runge     = zeros(NUM.NUMERICS.no_nodes_linear,2);
        PAR.dt_ini  = PAR.dt;
        
        for Runge = 1:1:4
            if Runge == 1 || Runge == 2
                PAR.dt = PAR.dt_ini/2;
            else
                PAR.dt = PAR.dt_ini;
            end
            [ NUM,MESH] = Solver(NUM,PAR,MESH,CHAR );
            
            Vx_runge(:,Runge) = NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes));
            Vz_runge(:,Runge) = NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes));
            P_runge(:,Runge)  = NUM.Solve.r(NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+1:end));
            
            MESH.GCOORD(1,:) = MESH.Old(1,:) + Vx_runge(:,Runge)'*PAR.dt;
            MESH.GCOORD(2,:) = MESH.Old(2,:) + Vz_runge(:,Runge)'*PAR.dt;
            
            display(['----- RK4 step ',num2str(Runge),' finished -----']) 
            
        end
        
        MESH = MESH_Old;
        
        PAR.dt = PAR.dt_ini;
        NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes))      = (1/6) * (Vx_runge(:,1)' + 2*Vx_runge(:,2)' + 2*Vx_runge(:,3)' + Vx_runge(:,4)');   % Vx update
        NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes))      = (1/6) * (Vz_runge(:,1)' + 2*Vz_runge(:,2)' + 2*Vz_runge(:,3)' + Vz_runge(:,4)');
        NUM.Solve.r(NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+1:end))  = (1/6) * (P_runge(:,1)' + 2*P_runge(:,2)' + 2*P_runge(:,3)' + P_runge(:,4)');
end

% Update strains for elasticity
[NUM] = Update_StressesAndStrains(NUM,PAR,MESH);




















