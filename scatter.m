
e0 = 8.854e-12;
speed_of_light = 299792458;
mu0 = 4 * pi * 1e-7;

frec = 300e6;
wavelength = speed_of_light / frec;
k0 = 2 * pi * frec / speed_of_light;
rot_frec = 2 * pi * frec;


% com_case, scat_case get values 1 to 3 to signify what radii to choose
comp_case = 3;
scat_case = 3;

switch scat_case
    case 1
        scatterer_radius = wavelength / 4;
    case 2
        scatterer_radius = wavelength / 2;
    case 3
        scatterer_radius = 5 * wavelength / 3;
    otherwise
        % assume case 2
        scatterer_radius = wavelength;
end
switch comp_case
    case 1
        comp_radius = wavelength / 2 + scatterer_radius;
    case 2
        comp_radius = wavelength + scatterer_radius;
    case 3
        comp_radius = 2 * wavelength + scatterer_radius;
end

% center coordinates
x0 = 0;
y0 = 0;

dx = 0.005;   % width of boundary regions

% Create the regions.....
in_bound_rad = scatterer_radius + dx;
out_bound_rad = comp_radius - dx;

%     comp. region   scatterer            inner boundary     outer boundary
gd = [1              1                    1                  1;
      x0             x0                   x0                 x0;
      y0             y0                   y0                 y0;
      comp_radius    scatterer_radius     in_bound_rad      out_bound_rad];

ns = char('comp','scatt','in_b','out_b')';
sf = 'comp - scatt';
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

Ei = zeros(Nn,1);
for id = 1:Nn
    Ei(id) = exp( -1i * k0 * p(1, id));
end


% node_id(i)=0 if node i has Dirichlet condition (known value), else it is 1.
node_id = ones(Nn,1);
Ez = zeros(Nn,1);
for id = 1:Nd
    if e(6,id) == 0 || e(7,id) == 0
        % one side of the edge is the outside of the mesh
        % the nodes are on the boundary
        if e(6,id) == 1 || e(7,id) == 1
            % nodes are in the inner boundary (region 1)
            % so they have known field = cos(k_0 x)
            node_id( e(1,id) ) = 0;
            node_id( e(2,id) ) = 0;
            Ez(e(1,id)) = - Ei(e(1,id));
            Ez(e(2,id)) = - Ei(e(2,id));
        end

    end
end


Nf = nnz(node_id);    % number of unknown nodes
Np = Nn - Nf;       % number of known nodes

% matrix of known values
Ez_p = zeros(Np,1);

% index(p) = new numbering of node p (position in Ez_f or Ez_p)
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
        Ez_p(counter_known) = Ez(ind);
    end
end


% calculate matrix Sff, Sfp, Tff, Tfp
Sff = spalloc(Nf,Nf,7*Nf);
Sfp = spalloc(Nf, Np, 7*Np);
Tff = spalloc(Nf,Nf,7*Nf);
Tfp = spalloc(Nf, Np, 7*Np);
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
            Sij = (1 / mu0) * (b(i)*b(j) + c(i)*c(j)) * Ae;
            if node_id(n(i)) == 1
                if node_id(n(j)) == 1
                    % both i and j are unknown
                    Sff(index(n(i)),index(n(j))) = Sff(index(n(i)),index(n(j))) + Sij;
                    if i == j
                        Tff(index(n(i)),index(n(j))) = Tff(index(n(i)),index(n(j))) + e0 * Ae / 6;
                    else
                        Tff(index(n(i)),index(n(j))) = Tff(index(n(i)),index(n(j))) + e0 * Ae / 12;
                    end
                else
                    % i is unknown and j is known
                    Sfp(index(n(i)), index(n(j))) = Sfp(index(n(i)), index(n(j))) + Sij;
                    if i == j
                        Tfp(index(n(i)),index(n(j))) = Tfp(index(n(i)),index(n(j))) + e0 * Ae / 6;
                    else
                        Tfp(index(n(i)),index(n(j))) = Tfp(index(n(i)),index(n(j))) + e0 * Ae / 12;
                    end
                end
            end
        end
    end
end

% find matrix C of the boundary matrix 
Cff = spalloc(Nf,Nf,7*Nf);
for edge = 1:Nd
    if e(6,edge) == 0 || e(7,edge) == 0
        % apply only on boundary edges
        node1 = e(1, edge);           
        node2 = e(2, edge);          
        x1 = p(1, node1);
        y1 = p(2, node1);
        x2 = p(1, node2);
        y2 = p(2, node2);
        length = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        Cff(index(node1),index(node2)) = Cff(index(node1),index(node2)) + length / (6 * mu0);
        Cff(index(node1),index(node1)) = Cff(index(node1),index(node1)) + length / (3 * mu0);
        Cff(index(node2),index(node2)) = Cff(index(node2),index(node2)) + length / (3 * mu0);
    end
end


alpha = -1i * k0 - 1 / (2 * comp_radius);

% solve system
Ez_f = (Sff - rot_frec .^ 2 * Tff - alpha * Cff) \ (- (Sfp - rot_frec .^ 2 * Tfp) * Ez_p );


% update Ez 
for ind = 1:Nn
    if node_id(ind) == 1  
        % update only for unknown nodes
        Ez(ind) = Ez_f(index(ind));
    end
end

field = Ei + Ez;

figure('Visible','off');
pdeplot(p,e,t,'XYData',abs(field)); axis equal tight;
colormap("jet");
exportgraphics(gcf, "./plots/scatter_case_"+scat_case+"_"+comp_case+".pdf", 'ContentType', 'vector');