

e0 = 8.854e-12;

% radii of inner conductor and dielectric
alpha = 0.76e-3;
beta = 1.75e-3;
% center coordinates
x0 = 0;
y0 = 0;

epsilon = 0.05e-3;   % width of boundary regions

%     conductor  dielectric         condactor bc       dielectric bc
gd = [1          1                  1                  1;
      x0        x0                 x0                 x0;
      y0        y0                 y0                 y0;
      alpha   beta   (alpha + epsilon)   (beta - epsilon)];


ns = char('cond','die','condBC','dieBC')';
sf = '(die-dieBC) + (dieBC-condBC) + (condBC-cond)';
dl = decsg(gd,sf,ns);



[p,e,t] = initmesh(dl);
[p,e,t] = refinemesh(dl,p,e,t);
[p,e,t] = refinemesh(dl,p,e,t);
%[p,e,t] = refinemesh(dl,p,e,t);


Nn = size(p,2);    % number of nodes
Ne = size(t,2);    % number of elements
Nd = size(e,2);    % number of edges


% node_id(i)=0 if node i has Dirichlet condition (known value), else it is 1.
node_id = ones(Nn,1);
% X0 contains for every node, the potential if it is known, else 0
X0 = zeros(Nn,1);
for id = 1:Nd
    if e(6,id) == 0 || e(7,id) == 0
    % one side of the edge is the outside of the mesh
    % the nodes are on the boundary
    node_id( e(1,id) ) = 0;
    node_id( e(2,id) ) = 0;

    if e(6,id) == 1 || e(7,id) == 1
        % nodes are in the inner boundary
        X0(e(1,id)) = 1;
        X0(e(2,id)) = 1;
    end

    end
end




Nf = nnz(node_id);    % number of unknown nodes
Np = Nn - Nf;       % number of known nodes

% matrix of known potentials
Fp = zeros(Np,1);

index = zeros(Nn,1);
counter_unknown = 0;
counter_known = 0;
for ind = 1:Nn
    if node_id(ind) == 1
        counter_unknown = counter_unknown + 1;
        index(ind) = counter_unknown;
    else
        counter_known = counter_known + 1;
        index(ind) = counter_known;
        Fp(counter_known) = X0(ind);
    end
end



% calculate matrix Sff and Sfp
Sff = spalloc(Nf,Nf,7*Nf);
Sfp = spalloc(Nf, Np, 7*Np);
for triangle = 1:Ne
    n(1:3) = t(1:3,triangle);
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3));
    D = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
    Ae = abs(D / 2);
    b(1) = (y(2)-y(3))/D;
    b(2) = (y(3)-y(1))/D;
    b(3) = (y(1)-y(2))/D;
    c(1) = (x(3)-x(2))/D;
    c(2) = (x(1)-x(3))/D;
    c(3) = (x(2)-x(1))/D;
    for i = 1:3
        for j = 1:3
            Sij = e0 * (b(i)*b(j) + c(i)*c(j)) * Ae;
            if node_id(n(i)) == 1
                if node_id(n(j)) == 1
                    % both i and j are unknown
                    Sff(index(n(i)),index(n(j))) = Sff(index(n(i)),index(n(j))) + Sij;
                else
                    % i is unknown and j is known
                    Sfp(index(n(i)), index(n(j))) = Sfp(index(n(i)), index(n(j))) + Sij;
                end
            end
        end
    end
end


Ff = Sff \ (-Sfp * Fp);


for ind = 1:Nn
    if node_id(ind) == 1
        X0(ind) = Ff(index(ind));
    end
end


% calculate the electric field E
[X0x, X0y] = pdegrad(p,t,X0);
Ex = -X0x;
Ey = -X0y;



% energy over unit length 
We = 0;
for triangle = 1:Ne
    n(1:3) = t(1:3,triangle);
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3));
    D = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
    Ae = abs(D / 2);
    b(1) = (y(2)-y(3))/D;
    b(2) = (y(3)-y(1))/D;
    b(3) = (y(1)-y(2))/D;
    c(1) = (x(3)-x(2))/D;
    c(2) = (x(1)-x(3))/D;
    c(3) = (x(2)-x(1))/D;
    for i = 1:3
        for j = 1:3
            We = We + 0.5 * e0 * X0(n(i)) * X0(n(j)) * (b(i)*b(j) + c(i)*c(j)) * Ae;
        end
    end
end

V = 1;
C =  2 * We / (V^2);

disp(We)
disp(C)



% plot regions
fig_reg = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15],"Visible","off");
pdegplot(dl, 'FaceLabels', 'on'); axis equal; axis tight;
title('regions');
exportgraphics(gcf, './plots/coaxial_regions.pdf', 'ContentType', 'vector');

% plot mesh
fig_mesh = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15],"Visible","off");
pdeplot(p,e,t); axis equal; axis tight;
title('triangulated mesh');
exportgraphics(gcf, './plots/coaxial_mesh.pdf', 'ContentType', 'vector');


% Plot the potential and field
fig_field = figure('Units', 'centimeters', 'Position', [1, 1, 15,15], 'Visible','off');
pdeplot(p,e,t,'XYData',X0,'FlowData',[Ex;Ey]); axis equal; axis tight;
title('potential and electric field')
colormap(winter);  
colorbar;
exportgraphics(gcf, './plots/coaxial_field.pdf', 'ContentType', 'vector');