%% 3D Stokes

clear all , close all
addpath('./AdditionalSource')


%% Constants

% Position parameters
W = 1;                
H = 1; 
D = 1;

rho_1 = 2;              
rho_2 = 1;

mu_1 = 1000;    
mu_2 = 100;
mu_3 = 10;

nphases = 11;  

h = [0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5];     

g = -10;                 
A0 = 0.05;            % maximum height of the sinus topo

% Numerical parameters
no_elems_x = 2;
no_elems_y = 2;
no_elems_z = 2;
no_elems_global = no_elems_x*no_elems_y*no_elems_z;

no_nodes_ele = 27;  
n_ele = no_nodes_ele;           % needed to reset no_nodes_ele for shape functions
no_nodes_ele_linear = 8;
no_intp = 27;       
ndof = 4;                       % 3 velocitys 1 pressure

no_nodes_x = (no_elems_x*2)+1;
no_nodes_y = (no_elems_y*2)+1;
no_nodes_z = (no_elems_z*2)+1;
no_nodes = no_nodes_x*no_nodes_y*no_nodes_z;

no_nodes_x_linear = (no_elems_x)+1;
no_nodes_y_linear = (no_elems_y)+1;
no_nodes_z_linear = (no_elems_z)+1;
no_nodes_linear = (no_elems_x+1)*(no_elems_y+1)*(no_elems_z+1);      


% Loop parameters
dx = (W/(no_nodes_x-1));
dy = (H/(no_nodes_y-1));
dz = (D/(no_nodes_z-1));

dt = 1;
nt = 30;

time = 0;
SecYear = 24*3600*365;

%% Get coordinates of the intp

[ weight,COORD ] = INT_PROPS(no_intp) ;     % get coordinates from the database with respect to the node/INTP system


%% initialize Model

% initialize x and z
x = [0:dx:W];
y = [0:dy:H];
z = [0:dz:D];
% x = x + 0.2*[rand(size(x))-0.5]*dx;             % Add random distribution additionally for variable spacing
% y = y + 0.2*[rand(size(y))-0.5]*dy;             % Add random distribution additionally for variable spacing
% z = z + 0.2*[rand(size(z))-0.5]*dz;             % Add random distribution additionally for variable spacing
[X2d,Y2d,Z2d] = meshgrid(x,y,z);

GCOORD = zeros(3,no_nodes_y*no_nodes_z);
i1 = 0;
i2 = 0;
for j = 1:1:no_nodes_z
    j1 = ((j-1)*no_nodes_x*no_nodes_y)+1;
    j2 = ((j)*no_nodes_x*no_nodes_y);
    GCOORD(3,j1:j2) = z(j);
    for i=1:1:no_nodes_y
        i1 = ((i-1)*no_nodes_x)+1 + ((j-1)*(no_nodes_x*no_nodes_y));
        i2 = ((i)*no_nodes_x)+ ((j-1)*(no_nodes_x*no_nodes_y));
        GCOORD(1,i1:i2) = x;
        GCOORD(2,i1:i2) = y(i);
    end

end

rho = zeros(size(GCOORD(2,:)));
mu  = zeros(size(GCOORD(2,:)));

% % Set phases and properties
% % Phase 1
% ind = find(GCOORD(2,:)<=H  & GCOORD(2,:)>=h(1,1));
% rho(ind) = rho_1;
% mu(ind) = mu_1;
% 
% % Phase 2
% ind = find(GCOORD(2,:)<h(1,1) & GCOORD(2,:)>=h(1,2));
% rho(ind) = rho_1;
% mu(ind) = mu_2;
% 
% % Phase 3
% ind = find(GCOORD(2,:)<h(1,2) & GCOORD(2,:)>=h(1,3));
% rho(ind) = rho_1;
% mu(ind) = mu_1;
% 
% % Phase 4
% ind = find(GCOORD(2,:)<h(1,3) & GCOORD(2,:)>=h(1,4));
% rho(ind) = rho_1;
% mu(ind) = mu_2;
% 
% % Phase 5
% ind = find(GCOORD(2,:)<h(1,4) & GCOORD(2,:)>=h(1,5));
% rho(ind) = rho_1;
% mu(ind) = mu_1;
% 
% % Phase 6
% ind = find(GCOORD(2,:)<h(1,5) & GCOORD(2,:)>=h(1,6));
% rho(ind) = rho_1;
% mu(ind) = mu_2;
% 
% % Phase 7
% ind = find(GCOORD(2,:)<h(1,6) & GCOORD(2,:)>=h(1,7));
% rho(ind) = rho_1;
% mu(ind) = mu_1;
% 
% % Phase 8
% ind = find(GCOORD(2,:)<h(1,7) & GCOORD(2,:)>=h(1,8));
% rho(ind) = rho_1;
% mu(ind) = mu_2;
% 
% % Phase 9
% ind = find(GCOORD(2,:)<h(1,8) & GCOORD(2,:)>=h(1,9));
% rho(ind) = rho_1;
% mu(ind) = mu_1;
% 
% % Phase 10
% ind = find(GCOORD(2,:)<h(1,9) & GCOORD(2,:)>=h(1,10));
% rho(ind) = rho_1;
% mu(ind) = mu_2;
% 
% % Phase 11
% ind = find(GCOORD(2,:)<h(1,10));
% rho(ind) = rho_2;
% mu(ind) = mu_3;

% Phase 11
ind = find(GCOORD(2,:)<h(1,10));
rho(ind) = rho_2;
mu(ind) = mu_2;
% Phase 11
ind = find(GCOORD(2,:)>=h(1,10));
rho(ind) = rho_1;
mu(ind) = mu_3;

% % make sinus in the middle of the model
% st                              = A0/sin(pi/2);
% i                               = (no_nodes_x*no_nodes_y):-no_nodes_x:no_nodes_x;
% factor                          = zeros(1,length(i));
% 
% factor(1,1:(length(i)-1)/2)     = linspace(0,st,(length(i)-1)/2);
% factor(1,(length(i)+1)/2:end)   = linspace(st,0,(length(i)+1)/2);
% 
% j = length(GCOORD)+1-no_nodes_x;
% i = length(GCOORD);
% 
% for l = 1:1:no_nodes_z
%     k = 1;
%     for o = (no_nodes_x*no_nodes_y):-no_nodes_x:no_nodes_x
%         k_topo          = linspace((0),(pi),no_nodes_x);
%         %     k_topo          = k_topo + 0.2*[rand(size(k_topo))-0.5];
%         z_topo          = [sin(k_topo)];
%         z_topo          = factor(1,k).*z_topo;
%         GCOORD(2,j:i)   = GCOORD(2,j:i) + [z_topo];
%         j               = j-no_nodes_x;
%         i               = i-no_nodes_x;
%         k               = k+1;
%     end
% end


% make random disribution for the middle of the model
% bounding = 2*H/10;
% ind = find(GCOORD(2,:)<h(1,10)+bounding & GCOORD(2,:)>h(1,10)-bounding);
% GCOORD(2,ind) = GCOORD(2,ind) + 0.06*[rand(size(GCOORD(2,ind)))];




%% Numbering
% UNIVESAL NUMBERINGS
% numbering degress of freedom for velocity
k = 0;
for i = 1:no_nodes
    for j = 1:ndof-1
        k = k+1;
        number_dof(j,i) = k;
    end
end


% numbering degrees of freedom for pressure
k = max(number_dof,[],2)+1;
k = k(2,1);
for i = no_nodes+1:no_nodes+1+no_nodes_linear-1
    for j = 1:ndof-3
        k = k+1;
        number_dof(j,i) = k;
    end
end


% numbering global nodes quadratic
j = 0;
for o = 1:1:no_nodes_z
    for n = no_nodes_y:-1:1
        for m = 1:1:no_nodes_x
            j = j+1;
            number_3d_quad(n,m,o) = j;
        end
    end
end

% numbering global nodes linear
j = 1;
i = 1;
k = 1;
for o = 1:2:no_nodes_z
    j = 1;
    for n = 1:2:no_nodes_y
        i = 1;
        for m = 1:2:no_nodes_x
            number_3d_linear(j,i,k) = number_3d_quad(n,m,o);
            i = i+1;
        end
        j = j+1;
    end
    k = k+1;
end

% numbering global for pressure
j = 0;
for o = 1:1:no_nodes_z_linear
    for n = no_nodes_y_linear:-1:1
        for m = 1:1:no_nodes_x_linear
            j = j+1;
            number_pressure(n,m,o) = j;
        end
    end
end

% numbering global degrees of freedom (for boundaries)
[n,m,o] = size(number_3d_linear);
for i=1:o
number_3d_linear_ud(:,:,i) = flipdim(number_3d_linear(:,:,i),i)';
end

q = number_dof(1,1:no_nodes);
w = number_dof(2,1:no_nodes);
r = number_dof(3,1:no_nodes);
e(number_3d_linear_ud(:)) = number_dof(1,no_nodes+1:no_nodes+no_nodes_linear);
j = 1;
k = 1;
for o = 1:1:no_nodes_z
    for n = no_nodes_y:-1:1
        for m = 1:1:no_nodes_x
            number_3d_ele_dof(n,m,o,1) = q(1,j);
            number_3d_ele_dof(n,m,o,2) = w(1,j);
            number_3d_ele_dof(n,m,o,3) = r(1,j);
            number_3d_ele_dof(n,m,o,4) = e(1,j);
            j = j+1;
        end
    end
end

% NON-UNIVERSAL NUMBERINGS
% numbering quadratic element
[k,l,p] = size(number_3d_quad);
j = 0;
u = no_nodes_x*no_nodes_y;          % third dimension operator to go to the next one you mist multiply k by 2 or 3
k = u;
for n = 1:1:no_elems_global
    if mod(n,(no_elems_x)) == 1 && n~=1
        j = j+l+1;
    end
    if mod(n,(no_elems_x*no_elems_y)) == 1 && n~=1
        k = k+2*u;
        j = 0;
        j = (j+(k-u))-n+1;
    end
    number_quad(1,n) = n+j;
    number_quad(2,n) = n+j+2;
    number_quad(3,n) = n+j+(l*2)+2;
    number_quad(4,n) = n+j+(l*2);
    number_quad(5,n) = 2*u+n+j;
    number_quad(6,n) = 2*u+n+j+2;
    number_quad(7,n) = 2*u+n+j+(l*2)+2;
    number_quad(8,n) = 2*u+n+j+(l*2);
    number_quad(9,n) = u+n+j+l+1;
    number_quad(10,n) = u+n+j;
    number_quad(11,n) = u+n+j+(l*2);
    number_quad(12,n) = n+j+l;
    number_quad(13,n) = 2*u+n+j+l;
    number_quad(14,n) = n+j+1;
    number_quad(15,n) = 2*u+n+j+1;
    number_quad(16,n) = n+j+(l*2)+1;
    number_quad(17,n) = 2*u+n+j+(l*2)+1;
    number_quad(18,n) = u+n+j+2;
    number_quad(19,n) = u+n+j+(l*2)+2;
    number_quad(20,n) = n+j+l+2;
    number_quad(21,n) = 2*u+n+j+l+2;
    number_quad(22,n) = u+n+j+l;
    number_quad(23,n) = u+n+j+l+2;
    number_quad(24,n) = u+n+j+1;
    number_quad(25,n) = u+n+j+(l*2)+1;
    number_quad(26,n) = n+j+l+1;
    number_quad(27,n) = 2*u+n+j+l+1;
    j = j+1;
end

% numbering for elemts and equations

[k,l_quad,p] = size(number_3d_quad);
j = 0;
o = 0;
u = no_nodes_x*no_nodes_y;          % third dimension operator to go to the next one you must multiply k by 1 or 2
v = no_nodes_x_linear*no_nodes_y_linear;
l_lin = no_nodes_x_linear;
k = u;       % index for quad 3D
h = v;       % index for linear 3D
for n = 1:1:no_elems_global
    if mod(n,(no_elems_x)) == 1 && n~=1
        j = j+l_quad+1;
        o = o+1;
    end
    if mod(n,(no_elems_x*no_elems_y)) == 1 && n~=1
        k = k+2*u;
        h = h+v;
        j = 0;                      % reset j and o
        j = (j+(k-u))-n+1;
        o = 0;
        o = (o+(h-v))-n+1;
    end
    number_ele_dof(1:3,n)   = number_dof(:,n+j);
    number_ele_dof(4:6,n)   = number_dof(:,n+j+2);
    number_ele_dof(7:9,n)   = number_dof(:,n+j+(l_quad*2)+2);
    number_ele_dof(10:12,n)   = number_dof(:,n+j+(l_quad*2));
    number_ele_dof(13:15,n)  = number_dof(:,2*u+n+j);
    number_ele_dof(16:18,n) = number_dof(:,2*u+n+j+2);
    number_ele_dof(19:21,n) = number_dof(:,2*u+n+j+(l_quad*2)+2);
    number_ele_dof(22:24,n) = number_dof(:,2*u+n+j+(l_quad*2));
    number_ele_dof(25:27,n) = number_dof(:,u+n+j+l_quad+1);
    number_ele_dof(28:30,n) = number_dof(:,u+n+j);
    number_ele_dof(31:33,n) = number_dof(:,u+n+j+(l_quad*2));
    number_ele_dof(34:36,n) = number_dof(:,n+j+l_quad);
    number_ele_dof(37:39,n) = number_dof(:,2*u+n+j+l_quad);
    number_ele_dof(40:42,n) = number_dof(:,n+j+1);
    number_ele_dof(43:45,n) = number_dof(:,2*u+n+j+1);
    number_ele_dof(46:48,n) = number_dof(:,n+j+(l_quad*2)+1);
    number_ele_dof(49:51,n) = number_dof(:,2*u+n+j+(l_quad*2)+1);
    number_ele_dof(52:54,n) = number_dof(:,u+n+j+2);
    number_ele_dof(55:57,n) = number_dof(:,u+n+j+(l_quad*2)+2);
    number_ele_dof(58:60,n) = number_dof(:,n+j+l_quad+2);
    number_ele_dof(61:63,n) = number_dof(:,2*u+n+j+l_quad+2);
    number_ele_dof(64:66,n) = number_dof(:,u+n+j+l_quad);
    number_ele_dof(67:69,n) = number_dof(:,u+n+j+l_quad+2);
    number_ele_dof(70:72,n) = number_dof(:,u+n+j+1);
    number_ele_dof(73:75,n) = number_dof(:,u+n+j+(l_quad*2)+1);
    number_ele_dof(76:78,n) = number_dof(:,n+j+l_quad+1);
    number_ele_dof(79:81,n) = number_dof(:,2*u+n+j+l_quad+1);
    number_ele_dof(82,n) = number_dof(1,no_nodes+n+o);
    number_ele_dof(83,n) = number_dof(1,no_nodes+n+o+1);
    number_ele_dof(84,n) = number_dof(1,no_nodes+n+o+l_lin+1);
    number_ele_dof(85,n) = number_dof(1,no_nodes+n+o+l_lin);
    number_ele_dof(86,n) = number_dof(1,no_nodes+n+o+v);
    number_ele_dof(87,n) = number_dof(1,no_nodes+n+o+v+1);
    number_ele_dof(88,n) = number_dof(1,no_nodes+n+o+v+l_lin+1);
    number_ele_dof(89,n) = number_dof(1,no_nodes+n+o+v+l_lin);
    j = j+1;
end






%% Boundaries / velocity field

% Constant pushing in x direction from the right and free surface
ux = W/100;
uy = 0;
uz = 0;
bcdof1  = [ number_3d_ele_dof(:,1,:,1)];     % left vx
bcdof2  = [number_3d_ele_dof(:,end,:,1)];    % right vx
bcdof3  = [number_3d_ele_dof(:,1,:,2)];      % left vy
bcdof4  = [number_3d_ele_dof(:,end,:,2)];    % right vy
bcdof5  = [number_3d_ele_dof(:,1,:,3)];      % left vz
bcdof6  = [number_3d_ele_dof(:,end,:,3)];    % right vz
bcdof7  = [number_3d_ele_dof(end,:,:,2)];    % bottom boundary vy
% bcdof8  = [number_3d_ele_dof(end,:,:,3)];    % bottom boundary vz
bcdof9  = [number_3d_ele_dof(:,:,1,3)];    % front boundary vy
% bcdof10 = [number_3d_ele_dof(:,:,1,3)];    % front boundary vz
bcdof11 = [number_3d_ele_dof(:,:,end,3)];    % back boundary vy
% bcdof12 = [number_3d_ele_dof(:,:,end,3)];    % back boundary vz


bcdof1  = bcdof1(:);
bcdof2  = bcdof2(:);
bcdof3  = bcdof3(:);
bcdof4  = bcdof4(:);
bcdof5  = bcdof5(:);
bcdof6  = bcdof6(:);
bcdof7  = bcdof7(:);
% bcdof8  = bcdof8(:);
bcdof9  = bcdof9(:);
% bcdof10 = bcdof10(:);
bcdof11 = bcdof11(:);
% bcdof12 = bcdof12(:);

% bcdof = [bcdof1;bcdof2;bcdof3;bcdof4;bcdof5;bcdof6;bcdof7;bcdof8;bcdof9;bcdof10;bcdof11;bcdof12];
bcdof = [bcdof1;bcdof2;bcdof3;bcdof4;bcdof5;bcdof6;bcdof7;bcdof9;bcdof11];

% bcval = [zeros(size(bcdof1)); ones(size(bcdof2))*(-ux); zeros(size(bcdof3)); ones(size(bcdof4))*(-uy); zeros(size(bcdof5)); ones(size(bcdof6))*(-uz); zeros(size(bcdof7)) ; zeros(size(bcdof8)); zeros(size(bcdof9)); zeros(size(bcdof10)); zeros(size(bcdof11)); zeros(size(bcdof12)) ];
bcval = [zeros(size(bcdof1)); ones(size(bcdof2))*(-ux); zeros(size(bcdof3)); ones(size(bcdof4))*(-uy); zeros(size(bcdof5)); ones(size(bcdof6))*(-uz); zeros(size(bcdof7)) ; zeros(size(bcdof9)); zeros(size(bcdof11)) ];
bcval = bcval(:);

m = [1;1;1;0;0;0];
 
%% initilize globals

% initialize global matrices

L = sparse(no_nodes*(ndof-1)+no_nodes_linear,no_nodes*(ndof-1)+no_nodes_linear);
FG = zeros(no_nodes*(ndof-1)+no_nodes_linear,1);
r = zeros(no_nodes*(ndof-1)+no_nodes_linear,1);
B  = zeros(6,no_nodes_ele*(ndof-1));
F  = zeros(no_nodes_ele*(ndof-1)+no_nodes_ele_linear,1);
D = eye(6,6);

%% Matrix building

for i = 1:no_elems_global
    
    % initialize local matrices
    KM = zeros(no_nodes_ele*(ndof-1),no_nodes_ele*(ndof-1));
    F  = zeros(no_nodes_ele*(ndof-1)+no_nodes_ele_linear,1);
    GM = zeros(no_nodes_ele*(ndof-1),no_nodes_ele_linear);
    LM = zeros(no_nodes_ele*(ndof-1)+no_nodes_ele_linear,no_nodes_ele*(ndof-1)+no_nodes_ele_linear);
    
    
    LGCOORD = (GCOORD(:,number_quad(:,i)))';      % local global coordinates for this element
    

    
    
    
    
    for j =1:no_intp
        
        [ NI,dNI ] = Shape_Functions( no_nodes_ele,j,COORD );

        
        N   = NI;
        dNds = dNI;
        J    = dNds*LGCOORD;
        
        % Inverse of Jacobian
        invJ = inv(J);
        detJ = det(J);
        dNdX = invJ*dNds;
        
        
        B(1,1:3:end-2) = dNdX(1,:);
        B(2,2:3:end-1)   = dNdX(2,:);
        B(3,3:3:end) = dNdX(3,:);
        B(4,1:3:end-2) = dNdX(2,:);
        B(4,2:3:end-1) = dNdX(1,:);
        B(5,2:3:end-1) = dNdX(3,:);
        B(5,3:3:end) = dNdX(2,:);
        B(6,1:3:end-2) = dNdX(3,:);
        B(6,3:3:end) = dNdX(1,:);
        
        D = [   2*mu(number_quad(9,i)) 0 0 0 0 0
                0 2*mu(number_quad(9,i)) 0 0 0 0
                0 0 2*mu(number_quad(9,i)) 0 0 0 
                0 0 0 mu(number_quad(9,i)) 0 0
                0 0 0 0 mu(number_quad(9,i)) 0
                0 0 0 0 0 mu(number_quad(9,i)) ];
        
        
        KM = KM + B'*D*B * weight(j)*detJ;
        F(2:3:end-8,1) = F(2:3:end-8,1) + (N'  .* (rho(number_quad(:,i))') .* g) * weight(j) * detJ;           % vy in F      % -8 because the last 8 values are zero due to the pressure
     
        % Compute linear shape functions
        no_nodes_ele = 8;
        [ NI,dNI ] = Shape_Functions( no_nodes_ele,j,COORD );
        N_p    = NI;
        GM = GM + (B'*m*N_p* weight(j)*detJ);
        no_nodes_ele = n_ele;
 
    end
    
    
    LM(1:81,1:81) = KM;
    LM(1:81,82:89) = GM;
    LM(82:89,1:81) = GM';
    
    % add element stiffness matrix to global stiffness matrix
    L(number_ele_dof(:,i),number_ele_dof(:,i)) = L(number_ele_dof(:,i),number_ele_dof(:,i)) + LM;
    FG(number_ele_dof(:,i),1) = FG(number_ele_dof(:,i),1) + F;
    
end

% Add boundaries to the system
for i = 1:1:length(bcdof)
    
    L(bcdof(i),:) = 0;
    L(bcdof(i),bcdof(i)) = 1;
    FG(bcdof(i)) = bcval(i);
    
end



%% Timeloop

for n = 1:dt:nt
    
    % Compute r
    r = L\FG;
    
    % Compute r as matrix for visualisation
    r = r';
    Ux_vec = r(number_dof(1,1:no_nodes));
    Uy_vec = r(number_dof(2,1:no_nodes));
    Uz_vec = r(number_dof(3,1:no_nodes));
    P = r(number_dof(1,no_nodes+1:end));
    
    Ux_vec_3d = Ux_vec(number_3d_quad);
    Uy_vec_3d = Uy_vec(number_3d_quad);
    Uz_vec_3d = Uz_vec(number_3d_quad);
    P_2d = P(number_pressure);
    
    X = GCOORD(1,:);
    X = X(number_3d_quad);
    Y = GCOORD(2,:);
    Y = Y(number_3d_quad);
    Z = GCOORD(3,:);
    Z = Z(number_3d_quad);
    X_ele = GCOORD(1,:);
    X_ele = X_ele(number_3d_linear);
    Y_ele = GCOORD(2,:);
    Y_ele = Y_ele(number_3d_linear);
    Z_ele = GCOORD(3,:);
    Z_ele = Z_ele(number_3d_linear);
    
    rho_3d = rho(:);
    rho_3d = rho_3d(number_3d_quad);
    mu_3d = mu(:);
    mu_3d = mu_3d(number_3d_quad);
    
    
    
%     % Plot
%     figure(1)
%     sliceomatic(Ux_vec_3d)
%     colorbar
%     title('VX')
%     drawnow
%     
%     
%     figure(2)
%     sliceomatic(X,Y,Z,Uy_vec_3d)
%     colorbar
%     title('VY')
%     drawnow
%
%     figure(3)
%     sliceomatic(X,Y,Z,Uz_vec_3d)
%     colorbar
%     title('Vz')
%     drawnow

% pressure needs to be fixed the plotting

%     figure(4)
%     slice(X_ele,Y_ele,Z_ele,P_2d)
%     colorbar
%     title('Pressure')
%     drawnow

figure(5)
quiver3(X,Y,Z,Ux_vec_3d,Uy_vec_3d,Uz_vec_3d)
colorbar
title('Velocity plot')
xlabel('x - Länge')
ylabel('y - Höhe')
drawnow

%     figure(6)
%     sliceomatic1p0(X,Y,Z,rho_3d)
%     colorbar
%     title('Density')
%     drawnow
%
%     figure(7)
%     sliceomatic1p0(X,Y,Z,mu_3d)
%     colorbar
%     title('Viscosity')
%     drawnow

Topo = Y-H;
ind = find(Topo<0);
Topo(ind) = 0;


fname = ['Output',num2str(n)];
save(fname,'GCOORD','rho','mu','rho_3d','mu_3d','X','Y','Z','P','time','Ux_vec_3d','Uz_vec_3d','Uy_vec_3d');
fname = ['OutputVisc',num2str(n),'.vtk'];
vtkwrite(fname,'structured_grid',X,Y,Z,'scalars','Viscosity',mu_3d)
fname = ['OutputDens',num2str(n),'.vtk'];
vtkwrite(fname,'structured_grid',X,Y,Z,'scalars','Density',rho_3d)
fname = ['OutputTopo',num2str(n),'.vtk'];
vtkwrite(fname,'structured_grid',X,Y,Z,'scalars','Topo',Topo)
fname = ['OutputVel',num2str(n),'.vtk'];
vtkwrite(fname,'structured_grid',X,Y,Z,'vectors','Velocity',Ux_vec_3d,Uy_vec_3d,Uz_vec_3d)


% update time and GCOORD
GCOORD(1,:) = GCOORD(1,:) + Ux_vec*dt;
GCOORD(2,:) = GCOORD(2,:) + Uy_vec*dt;
GCOORD(3,:) = GCOORD(3,:) + Uz_vec*dt;

time = time+dt;
r = r';

    
end





