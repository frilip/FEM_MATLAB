%{

                  FEM method for capacitor (parallel)
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


We solve for a parallel capacitor of infinite length and conductors of 
width w, height h, distance between them d and dielectric between them
er.
The upper conductor has voltage +V/2 and the lower -V/2
The domain is square, relative to w and can be modified.

The mesh is created using [p,e,t] = initmesh()

The global numbering of nodes is given by their position in the p matrix
(matrix of all nodes)
the local by their position in the t matrix (triangle)

for more, read: 
https://www.mathworks.com/help/pde/ug/mesh-data-pet-triples.html


We will also need a map from the numbering given in the matrix p 
to a different numbering where known and unknown nodes are differenciated
This is done by the index array


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

% dimentions
w = 0.04;
h = 0.002;
d = 0.01;
V = 100;
er = 2.2;

% computational domain is of dimention A x A 
domain_size_relative_to_w = 9;
A = domain_size_relative_to_w * w;

% create regions
frame = [3, 4, -A/2, A/2, A/2, -A/2, -A/2, -A/2, A/2, A/2];
conductor_down = [3, 4, -w/2, w/2, w/2, -w/2, -d/2-h, -d/2-h, -d/2, -d/2];
conductor_up = [3, 4, -w/2, w/2, w/2, -w/2, d/2, d/2, d/2+h, d/2+h];
dielectric = [3, 4, -w/2, w/2, w/2, -w/2, d/2, d/2, -d/2, -d/2];

% a region that surrounds the capacitor, will be used to identify dirichlet conditions 
surrounding_space = [3, 4, -w, w, w, -w, -d-2*h, -d-2*h, d+2*h, d+2*h];  


gd = [frame', conductor_down', conductor_up', dielectric',surrounding_space'];
 
ns = char('frame','cond_up','cond_down', 'dielectric','surr_sp')';
sf = 'frame - cond_up - cond_down';
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

    % if it is inside region 1 (surrounding the capacitor) 
    % or region 3 (in the cacacitor) then we have dirichlet conditions
        if e(6,id) == 1 || e(7,id) == 1 || e(6,id) == 3 || e(7,id) == 3
            % nodes belong to the capacitors edges
            node_id( e(1,id) ) = 0;
            node_id( e(2,id) ) = 0;
            
            if p(2,e(1,id)) > 0
                % y value above 0, belongs to the upper conductor (+V/2)
                % both of the nodes of the edge are in the same conductor
                % as the edge is on the boundary
                X0(e(1,id)) = V/2;
                X0(e(2,id)) = V/2;
            else
                % belongs to the lower conductor (-V/2)
                X0(e(1,id)) = -V/2;
                X0(e(2,id)) = -V/2;
            end
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
    region = t(4,triangle);
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
            % check if the triangle is inside the capacitor
            if region == 3
                % inside, use er * e0
                Sij = er * e0 * (b(i)*b(j) + c(i)*c(j)) * Ae;
            else
                % outside, assume air, use only e0
                Sij = e0 * (b(i)*b(j) + c(i)*c(j)) * Ae;
            end
            
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


% solve system
Ff = Sff \ (-Sfp * Fp);



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
    region = t(4,triangle);
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
            if region == 3
                % inside, use er * e0
                We = We + 0.5 * er * e0 * X0(n(i)) * X0(n(j)) * (b(i)*b(j) + c(i)*c(j)) * Ae;
            else
                % outside, assume air, use only e0
                We = We + 0.5 * e0 * X0(n(i)) * X0(n(j)) * (b(i)*b(j) + c(i)*c(j)) * Ae;
            end
        end
    end
end

% capacitance
C =  2 * We / (V^2);

fprintf('Using domain size %dw x %dw.',domain_size_relative_to_w, domain_size_relative_to_w);
fprintf('Used %d refinements. The degree of freedom (nodes with uknown potential) is: %d.\n', refine_amount, Nf);
fprintf('The calculated capacitance over unit length is: %d Farad/m.\n', C);




% Suppress warnings, because saving as vector type is slow, 
% and a warning saying that will appear
warning('off', 'all');


% The figures are not shown, but saved!!

% plot regions
fig_reg = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15], 'Visible','off');
pdegplot(dl, 'FaceLabels', 'on'); axis equal; axis tight;
title('regions');
exportgraphics(gcf, "./plots/capacitor_regions_domain:" + domain_size_relative_to_w + ".pdf", 'ContentType', 'vector');

% plot mesh
fig_mesh = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15], 'Visible','off');
pdeplot(p,e,t); axis equal; axis tight;
title("triangulated mesh, " + refine_amount + " refinements");
exportgraphics(gcf, "./plots/capacitor_mesh_" + refine_amount + "_domain:" + domain_size_relative_to_w + ".pdf", 'ContentType', 'vector');



% plot the potential
fig_potential = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15], 'Visible','off');
pdeplot(p,e,t,'XYData',X0); axis equal; axis tight;
title("potential, " + refine_amount + " refinements")
colormap(jet);  
colorbar;         
exportgraphics(gcf, "./plots/capacitor_potential_" + refine_amount + "_domain:" + domain_size_relative_to_w + ".pdf", 'ContentType', 'vector');



% Plot the field
fig_field = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15], 'Visible','off');
pdegplot(dl);
hold on;
pdeplot(p,e,t,'FlowData',[Ex;Ey]); axis equal; axis tight;
title("electric field, " + refine_amount + " refinements")      
exportgraphics(gcf, "./plots/capacitor_field_" + refine_amount + "_domain:" + domain_size_relative_to_w + ".pdf", 'ContentType', 'vector');





% Bonus, plot the field with streamlines
% Compute triangle centers
triCenters = (p(:, t(1,:)) + p(:, t(2,:)) + p(:, t(3,:))) / 3;

% Interpolate Ex and Ey (defined per triangle) onto a regular grid
F_Ex = scatteredInterpolant(triCenters(1,:)', triCenters(2,:)', Ex(:), 'linear', 'none');
F_Ey = scatteredInterpolant(triCenters(1,:)', triCenters(2,:)', Ey(:), 'linear', 'none');

% Create a uniform grid over the domain
[xg, yg] = meshgrid( ...
    linspace(min(p(1,:)), max(p(1,:)), 200), ...
    linspace(min(p(2,:)), max(p(2,:)), 200));

% Evaluate interpolated field
Ex_grid = F_Ex(xg, yg);
Ey_grid = F_Ey(xg, yg);

% Plot potential
fig_streamlines = figure('Units', 'centimeters', 'Position', [1, 1, 15, 15], 'Visible','off');
pdeplot(p, e, t, 'XYData', X0, 'Mesh', 'off'); 
hold on; axis equal tight;

% Plot electric field as streamlines (built-in, elegant)
h = streamslice(xg, yg, Ex_grid, Ey_grid, 2);  % '2' sets spacing
set(h, 'Color', 'k', 'LineWidth', 0.8);        % Make them black

% Aesthetics
title("Electric Potential and Field Lines, " + refine_amount + " refinements");
colormap(jet);
colorbar;

% Export
exportgraphics(gcf, "./plots/capacitor_field_streamlines_" + refine_amount + "_domain:" + domain_size_relative_to_w + ".pdf", 'ContentType', 'vector');
