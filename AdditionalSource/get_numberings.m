function [ NUM ] = get_numberings(NUM)


if NUM.NUMERICS.Dimension == 2
    %% Numbering
    
    % UNIVESAL NUMBERINGS
    % numbering degress of freedom for velocity
    k = 0;
    for i = 1:NUM.NUMERICS.no_nodes
        for j = 1:NUM.NUMERICS.ndof-1
            k = k+1;
            NUM.Number.number_dof(j,i) = k;
        end
    end
    NUM.Number.number_dof_vel = NUM.Number.number_dof;
    
    % extend numbering degrees of freedom for pressure
    NUM.NUMERICS.no_nodes_linear = (NUM.NUMERICS.no_nodes_x_linear)*(NUM.NUMERICS.no_nodes_z_linear); % attention midpoint is also linear
    k = max(NUM.Number.number_dof,[],2);
    k = k(2,1);
    for i = NUM.NUMERICS.no_nodes+1:NUM.NUMERICS.no_nodes+1+NUM.NUMERICS.no_nodes_linear-1
        for j = 1:NUM.NUMERICS.ndof-2
            k = k+1;
            NUM.Number.number_dof(j,i) = k;
        end
    end
    
    
    % numbering global nodes quadratic
    j = 0;
    for n = NUM.NUMERICS.no_nodes_z:-1:1
        for m = 1:1:NUM.NUMERICS.no_nodes_x
            j = j+1;
            NUM.Number.number_2d(n,m) = j;
        end
    end
    
    % numbering global nodes linear
    j = 1;
    i = 1;
    for n = 1:2:NUM.NUMERICS.no_nodes_z
        i = 1;
        for m = 1:2:NUM.NUMERICS.no_nodes_x
            NUM.Number.number_2d_linear(j,i) = NUM.Number.number_2d(n,m);
            i = i+1;
        end
        j = j+1;
    end
    
    % numbering global
    j = 0;
    for n = NUM.NUMERICS.no_elems_z:-1:1
        for m = 1:1:NUM.NUMERICS.no_elems_x
            j = j+1;
            NUM.Number.number_elements(n,m) = j;
        end
    end
    
    % numbering linear nodes
    j = 0;
    for n = NUM.NUMERICS.no_nodes_z_linear:-1:1
        for m = 1:1:NUM.NUMERICS.no_nodes_x_linear
            j = j+1;
            NUM.Number.number_pressure(n,m) = j;
        end
    end
    
    q = zeros(size(NUM.Number.number_dof(1,:)));
    w = zeros(size(NUM.Number.number_dof(1,:)));
    e = zeros(size(NUM.Number.number_dof(1,:)));
    
    q = NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes);
    w = NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes);
    e(1,NUM.NUMERICS.no_nodes+1:NUM.NUMERICS.no_nodes+NUM.NUMERICS.no_nodes_linear) = NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+1:NUM.NUMERICS.no_nodes+NUM.NUMERICS.no_nodes_linear);
    j = 1;
    k = 1;
    for n = NUM.NUMERICS.no_nodes_z:-1:1
        for m = 1:1:NUM.NUMERICS.no_nodes_x
            NUM.Number.number_2d_ele_dof(n,m,1) = q(1,j);
            NUM.Number.number_2d_ele_dof(n,m,2) = w(1,j);
            NUM.Number.number_2d_ele_dof(n,m,3) = e(1,j);
            j = j+1;
        end
    end
    
    % NON-UNIVERSAL NUMBERINGS
    switch NUM.NUMERICS.Elementtype
        case 'linear'
            % numbering linear element
            l = NUM.NUMERICS.no_nodes_x;
            j = 0;
            for n = 1:1:NUM.NUMERICS.no_elems_global
                if mod(n,(l-1)) == 1 && n~=1
                    j = j+1;
                end
                NUM.Number.number_quad(1,n) = n+j;
                NUM.Number.number_quad(2,n) = n+j+l;
                NUM.Number.number_quad(3,n) = n+j+l+1;
                NUM.Number.number_quad(4,n) = n+j+1;
            end
            
            % numbering for elemts and equations
            j = 0;
            l = NUM.NUMERICS.no_nodes_x;
            for n = 1:1:NUM.NUMERICS.no_elems_global
                if mod(n,(l-1)) == 1 && n~=1
                    j = j+1;
                end
                NUM.Number.number_ele_dof(1:2,n) = NUM.Number.number_dof(:,j + n);
                NUM.Number.number_ele_dof(3:4,n) = NUM.Number.number_dof(:,j+l+ n);
                NUM.Number.number_ele_dof(5:6,n) = NUM.Number.number_dof(:,j+l+1+ n);
                NUM.Number.number_ele_dof(7:8,n) = NUM.Number.number_dof(:,j+1+ n);
                NUM.Number.number_ele_dof(9,n)   = NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+n);
            end
            
        case 'quadratic'
            % numbering quadratic element
            l = NUM.NUMERICS.no_nodes_x;
            j = 0;
            for n = 1:1:NUM.NUMERICS.no_elems_global
                if mod(n,(NUM.NUMERICS.no_elems_x)) == 1 && n~=1
                    j = j+l+1;
                end
                NUM.Number.number_quad(1,n) = n+j;
                NUM.Number.number_quad(2,n) = n+j+l;
                NUM.Number.number_quad(3,n) = n+j+(l*2);
                NUM.Number.number_quad(4,n) = n+j+(l*2)+1;
                NUM.Number.number_quad(5,n) = n+j+(l*2)+2;
                NUM.Number.number_quad(6,n) = n+j+l+2;
                NUM.Number.number_quad(7,n) = n+j+2;
                NUM.Number.number_quad(8,n) = n+j+1;
                NUM.Number.number_quad(9,n) = n+j+l+1;
                j = j+1;
            end
            
            % numbering for elemts and equations
            j = 0;
            o = 0;
            l = NUM.NUMERICS.no_nodes_x;
            h = NUM.NUMERICS.no_elems_x+1;
            for n = 1:1:NUM.NUMERICS.no_elems_global
                if mod(n,(NUM.NUMERICS.no_elems_x)) == 1 && n~=1
                    j = j+l+1;
                    o = o+1;
                end
                NUM.Number.number_ele_dof(1:2,n)   = NUM.Number.number_dof(:,n+j);
                NUM.Number.number_ele_dof(3:4,n)   = NUM.Number.number_dof(:,n+j+l);
                NUM.Number.number_ele_dof(5:6,n)   = NUM.Number.number_dof(:,n+j+(l*2));
                NUM.Number.number_ele_dof(7:8,n)   = NUM.Number.number_dof(:,n+j+(l*2)+1);
                NUM.Number.number_ele_dof(9:10,n)  = NUM.Number.number_dof(:,n+j+(l*2)+2);
                NUM.Number.number_ele_dof(11:12,n) = NUM.Number.number_dof(:,n+j+l+2);
                NUM.Number.number_ele_dof(13:14,n) = NUM.Number.number_dof(:,n+j+2);
                NUM.Number.number_ele_dof(15:16,n) = NUM.Number.number_dof(:,n+j+1);
                NUM.Number.number_ele_dof(17:18,n) = NUM.Number.number_dof(:,n+j+l+1);
                NUM.Number.number_ele_dof(19,n) = NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+n+o);
                NUM.Number.number_ele_dof(20,n) = NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+n+o+h);
                NUM.Number.number_ele_dof(21,n) = NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+n+o+h+1);
                NUM.Number.number_ele_dof(22,n) = NUM.Number.number_dof(1,NUM.NUMERICS.no_nodes+n+o+1);
                j = j+1;
            end
    end
end

if NUM.NUMERICS.Dimension == 3
    error('function needs to be extended!')
end



end

