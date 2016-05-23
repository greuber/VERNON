function [ NUM ] = get_boundaries( NUM,boundary )

switch boundary
    case 'FreeSlipEverywhere'
        % Free slip everywhere
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)',  NUM.Number.number_2d_ele_dof(:,end,1)',  NUM.Number.number_2d_ele_dof(1,:,2),  NUM.Number.number_2d_ele_dof(end,:,2)];
        NUM.Boundary.bcval = [zeros(size(NUM.Number.number_2d_ele_dof(:,1,1)')), zeros(size(NUM.Number.number_2d_ele_dof(:,end,1)')), zeros(size(NUM.Number.number_2d_ele_dof(1,:,1))), zeros(size(NUM.Number.number_2d_ele_dof(end,:,1))) ];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
    case 'Rayleighbenchmark'
        % Free slip at the sides and no slip at bottom and top
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)',  NUM.Number.number_2d_ele_dof(:,end,1)',  NUM.Number.number_2d_ele_dof(1,:,2),NUM.Number.number_2d_ele_dof(1,:,1),  NUM.Number.number_2d_ele_dof(end,:,2),  NUM.Number.number_2d_ele_dof(end,:,1)];
        NUM.Boundary.bcval = [zeros(size(NUM.Number.number_2d_ele_dof(:,1,1)')), zeros(size(NUM.Number.number_2d_ele_dof(:,end,1)')), zeros(size(NUM.Number.number_2d_ele_dof(1,:,1))), zeros(size(NUM.Number.number_2d_ele_dof(1,:,1))), zeros(size(NUM.Number.number_2d_ele_dof(end,:,1))) , zeros(size(NUM.Number.number_2d_ele_dof(end,:,1)))];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
    case 'FreeSurface'
        % free surface at the top
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)', NUM.Number.number_2d_ele_dof(:,end,1)',NUM.Number.number_2d_ele_dof(end,:,2) ];
        NUM.Boundary.bcval = [zeros(size(NUM.Number.number_2d_ele_dof(:,1,1)')), zeros(size(NUM.Number.number_2d_ele_dof(:,end,1)')),zeros(size(NUM.Number.number_2d_ele_dof(end,:,1)))];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
    case 'ConstantPushing'
        % Constant pushing in x direction from the right and free slip at
        % the surface
        if isfield(NUM,'ebg')
            ux = NUM.ebg;
        elseif isfield(NUM,'vbg')
            ux = NUM.vbg;
        end
        uy = 0;
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)', NUM.Number.number_2d_ele_dof(:,end,1)', NUM.Number.number_2d_ele_dof(end,:,2),NUM.Number.number_2d_ele_dof(1,:,2)];
        NUM.Boundary.bcval = [zeros(size(NUM.Number.number_2d_ele_dof(:,1,1)')), ones(size(NUM.Number.number_2d_ele_dof(:,end,1)'))*(ux), zeros(size(NUM.Number.number_2d_ele_dof(:,1,2)')), ones(size(NUM.Number.number_2d_ele_dof(:,end,2)'))*(-uy), zeros(size(NUM.Number.number_2d_ele_dof(end,:,2))) ];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
    case 'ConstantPushing+'
        % Constant pushing in x direction from the right and free surface
        if isfield(NUM,'ebg')
            ux = NUM.ebg;
        elseif isfield(NUM,'vbg')
            ux = NUM.vbg;
        end
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)', NUM.Number.number_2d_ele_dof(:,end,1)', NUM.Number.number_2d_ele_dof(end,:,2)];
        NUM.Boundary.bcval = [zeros(size(NUM.Number.number_2d_ele_dof(:,1,1)')), ones(size(NUM.Number.number_2d_ele_dof(:,end,1)'))*(ux), zeros(size(NUM.Number.number_2d_ele_dof(end,:,2))) ];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
    case 'ConstantPushingboth'
        % free slip on the bottom and constant pushing from both sides
        if isfield(NUM,'ebg')
            ux = NUM.ebg;
        elseif isfield(NUM,'vbg')
            ux = NUM.vbg;
        end
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)', NUM.Number.number_2d_ele_dof(:,end,1)',  NUM.Number.number_2d_ele_dof(end,:,2)];
        NUM.Boundary.bcval = [ones(size(NUM.Number.number_2d_ele_dof(:,1,1)'))*(-ux), ones(size(NUM.Number.number_2d_ele_dof(:,end,1)'))*(ux), zeros(size(NUM.Number.number_2d_ele_dof(end,:,2)))];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
    case 'FoldingBenchmark'
        % Constant pushing in x direction from the right, free slip
        % elsewhere, free surface and no slip at the bottom
        if isfield(NUM,'ebg')
            ux = NUM.ebg;
        elseif isfield(NUM,'vbg')
            ux = NUM.vbg;
        end
        NUM.Boundary.bcdof = [ NUM.Number.number_2d_ele_dof(:,1,1)', NUM.Number.number_2d_ele_dof(:,end,1)', NUM.Number.number_2d_ele_dof(end,:,2), NUM.Number.number_2d_ele_dof(end,:,1) ];
        NUM.Boundary.bcval = [zeros(size(NUM.Number.number_2d_ele_dof(:,1,1)')), ones(size(NUM.Number.number_2d_ele_dof(:,end,1)'))*(ux), zeros(size(NUM.Number.number_2d_ele_dof(end,:,2))), zeros(size(NUM.Number.number_2d_ele_dof(end,:,1))) ];
        NUM.Boundary.bcdof = NUM.Boundary.bcdof(:);
        NUM.Boundary.bcval = NUM.Boundary.bcval(:);
        
end

