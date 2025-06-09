%{

                  FEM method for coaxial cable
*-------------------------------------------------------------------------*

In the FEM method we triangulate the mesh and estimate the unknown function
on the nodes of the triangles. These values are then interpolated
to find the value anywhere on the mesh.

The FEM method for solving potential for 2d electrostatic problems with
Dirichlet or Neumann boundary conditions equates to solving the linear 
system:
                    Sff * Ff = - Sfp * Fp 
where:
* Ff is the column vector that contains all unknown potentials (on nodes)
* Fp is the column vector that contains all known potentials
and S is calculated by

if for nodes i,j (local numbering) in a triangle, with coordinates 
(xi,yi), (xj,yj), i,j in {1,2,3}:
- Ae is the area of the triangle 
- D = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
- b1 = (y(2)-y(3))/D;
- c1 = (x(3)-x(2))/D;
- bi,ci are given with circular rotation from above equations
then:

for nodes p,q (global numbering) with local numberin i,j 
S(p,q) = sum e * (bibj + cicj) * Ae, for triangles where both i and j
belong
(e is the dielectric constant) 

* Sff = S(p,q) for both p and q unknown 
* Sfp = S(p,q) for p unknown and q known

Sff and Sfp are 0 if p and q are not neighboars, so they are sparse
matrixes.
They can easily be calculated by:

For every triangle:
    For every combination of nodes i,j (local) -> p,q (global):
        S(p,q) += e * (bibj + cicj) * Ae


We solve for a coaxial cable with inner and outer radii: a and b.
The mesh is created using [p,e,t] = initmesh()

The global numbering of nodes is given by their position in the p matrix
(matrix of all nodes)
the local by their position in the t matrix (triangle)

for more, read: 
https://www.mathworks.com/help/pde/ug/mesh-data-pet-triples.html


We will also need a map from the numbering given in the matrix p 
to a different numbering where known and unknown nodes are differenciated
This is done by the index array

Direct and iterative solvers will be used. Benchmarks will be printed.


The energy over unit length will be calculated by adding:
      0.5 * e0 * X0(n(i)) * X0(n(j)) * (b(i)*b(j) + c(i)*c(j)) * Ae;
where X0(n(i)) is the potential (estimated) of node i, for every
combination of nodes in every triangle

for more info and derivation of the above see the report (in same
directory)

The capacitanve over unit length is calculated from the energy.


The graphs are saved in vector form, which makes the procedure slow
(some seconds),
if you need faster results, remove 'ContentType', 'vector' 
from exportgraphics.
%}

e0 = 8.854e-12;

% radii of inner conductor and dielectric
alpha = 0.76e-3;
beta = 1.75e-3;
% center coordinates
x0 = 0;
y0 = 0;

dx = 0.05e-3;   % width of boundary regions

% actual capacitance (calculated analytically) is: 
analytic_capacitance = 6.67014293e-11;  

% Create the regions.....

%     conductor  dielectric         condactor bc       dielectric bc
gd = [1          1                  1                  1;
      x0        x0                 x0                 x0;
      y0        y0                 y0                 y0;
      alpha   beta       (alpha + dx)         (beta - dx)];


ns = char('cond','die','condBC','dieBC')';
sf = '(die-dieBC) + (dieBC-condBC) + (condBC-cond)';
dl = decsg(gd,sf,ns);


% create triangular mesh
[p,e,t] = initmesh(dl);
% refine 
refine_amount = 2;
for i = 1:refine_amount
    [p,e,t] = refinemesh(dl,p,e,t);
end



Nn = size(p,2);    % number of nodes
Ne = size(t,2);    % number of elements
Nd = size(e,2);    % number of edges


% initialise array node_id and X0 
% both Nn x 1, node_id helps us tell which nodes have known and unkown
% values
% X0 is the potential array, it is initialised as 0 everywhere except
% where the potential is known (Dirichlet boundary conditions)

% node_id(i)=0 if node i has Dirichlet condition (known value), else it is 1.
node_id = ones(Nn,1);
X0 = zeros(Nn,1);
for id = 1:Nd
    if e(6,id) == 0 || e(7,id) == 0
        % one side of the edge is the outside of the mesh
        % the nodes are on the boundary
        node_id( e(1,id) ) = 0;
        node_id( e(2,id) ) = 0;
    
        if e(6,id) == 1 || e(7,id) == 1
            % nodes are in the inner boundary (region 1)
            % so they have known potential 1 Volt
            X0(e(1,id)) = 1;
            X0(e(2,id)) = 1;
        end

    end
end




Nf = nnz(node_id);    % number of unknown nodes
Np = Nn - Nf;       % number of known nodes

% matrix of known potentials
Fp = zeros(Np,1);


% index(p) = new numbering of node p (position in Ff or Fp)
% the code below counts the number of known and unknown nodes encountered
% set index to the corresponding value

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


%---------------------%
% Direct Solver
%---------------------%
tic;
Ff_direct = Sff \ (-Sfp * Fp);
time_direct = toc;
fprintf('Direct solver time: %.6f seconds\n', time_direct);

%---------------------%
% Iterative Solver 
%---------------------%
tic;
tol = 0.009;
max_it = 150;
Ff_pcg = pcg(Sff, (-Sfp * Fp), tol, max_it);
time_pcg = toc;
fprintf ('pcg solver time: %.6f seconds\n',time_pcg);



% choose solver 
solver = "Direct";  % change this variable to chose different solver

if solver == "pcg"
    Ff =  Ff_pcg;
    fprintf('Using iterative solver!!\n');
elseif solver == "Direct"
    Ff = Ff_direct;
    fprintf('Using direct solver!!\n')
else
    disp("invalid solver, going with direct...");
    Ff = Ff_direct;
end


% update X0 to hold calculated potential 
for ind = 1:Nn
    if node_id(ind) == 1  
        % update only for unknown nodes
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


% capacitance
V = 1;
C =  2 * We / (V^2);

fprintf('Used %d refinements. The degree of freedom (nodes with uknown potential) is: %d.\n', refine_amount, Nf);
fprintf('The calculated capacitance over unit length is: %d Farad/m.\n', C);
fprintf('capacitance relative error is: %d%%.\n', abs((C - analytic_capacitance) / analytic_capacitance) * 100 );


% Suppress warnings, because saving as vector type is slow, 
% and a warning saying that will appear
warning('off', 'all');

% The figures are not shown, but saved!!

% plot regions
fig_reg = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15],"Visible","off");
pdegplot(dl, 'FaceLabels', 'on'); axis equal; axis tight;
title("regions");
exportgraphics(gcf, './plots/coaxial_regions.pdf', 'ContentType', 'vector');

% plot mesh
fig_mesh = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15],"Visible","off");
pdeplot(p,e,t); axis equal; axis tight;
title("triangulated mesh, " + refine_amount + " refinements");
exportgraphics(gcf, "./plots/coaxial_mesh_"+refine_amount+".pdf", 'ContentType', 'vector');


% Plot the potential and field
fig_field = figure('Units', 'centimeters', 'Position', [1, 1, 15,15], 'Visible','off');
pdeplot(p,e,t,'XYData',X0,'FlowData',[Ex;Ey]); axis equal; axis tight;
title("      potential and electric field, " + refine_amount + " refinements, " + solver + " solver");
colormap(cool);  
colorbar;
exportgraphics(gcf, "./plots/coaxial_field_"+refine_amount+"_"+solver+".pdf", 'ContentType', 'vector');