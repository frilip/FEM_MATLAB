function plot_streamlines(Fx, Fy, p, e, t, color)
% plot_streamlines_only(Fx, Fy, p, e, t, color)
% Plots only streamlines of the vector field F = (Fx, Fy),
% where Fx and Fy are defined per triangle.
%
% Parameters:
%   Fx, Fy - vector field components (per triangle)
%   p, e, t - mesh data
%   color - streamline color (e.g., 'k', 'r', [0 0.5 1])

    if nargin < 6
        color = 'k';  % Default color is black
    end

    % Compute triangle centers
    triCenters = (p(:, t(1,:)) + p(:, t(2,:)) + p(:, t(3,:))) / 3;

    % Interpolate field onto grid
    F_Ex = scatteredInterpolant(triCenters(1,:)', triCenters(2,:)', Fx(:), 'linear', 'none');
    F_Ey = scatteredInterpolant(triCenters(1,:)', triCenters(2,:)', Fy(:), 'linear', 'none');

    % Create grid
    [xg, yg] = meshgrid( ...
        linspace(min(p(1,:)), max(p(1,:)), 200), ...
        linspace(min(p(2,:)), max(p(2,:)), 200));

    Ex_grid = F_Ex(xg, yg);
    Ey_grid = F_Ey(xg, yg);

    % Plot streamlines only
    h = streamslice(xg, yg, Ex_grid, Ey_grid, 1); % Adjust spacing for density
    set(h, 'Color', color, 'LineWidth', 0.4);

    axis equal tight;
end






% ------------------- MAIN ---------------------------





e0 = 8.854e-12;
mu0 = 4 * pi * 1e-7;
speed_of_light = 299792458;

radius = 0.01;
% center coordinates
x0 = 0;
y0 = 0;

% create the regions..

gd = [1;
      x0;
      y0;
      radius];

ns = char('waveguide')';
sf = 'waveguide';
dl = decsg(gd, sf, ns);

% create triangular mesh
[p,e,t] = initmesh(dl);
% refine 
refine_amount = 3;
for i = 1:refine_amount
    [p,e,t] = refinemesh(dl,p,e,t);
end

Nn = size(p,2);    % number of nodes
Ne = size(t,2);    % number of elements
Nd = size(e,2);    % number of edges

%{
TE:
    Solve for Hz, with homogenous Neumann boundary conditions
TM:
    Solve for Ez, with Dirichlet Ez=0 boundary conditions

So TE doesn't have known nodes, and TM does.
%}

Ez_TM = zeros(Nn,1);
Hz_TE = zeros(Nn,1);

% find known nodes of Ez
% node_id(i)=0 if node i has Dirichlet condition (known value), else it is 1.
node_id = ones(Nn,1);
for id = 1:Nd
    if e(6,id) == 0 || e(7,id) == 0
        % one side of the edge is the outside of the mesh
        % the nodes are on the boundary
        node_id( e(1,id) ) = 0;
        node_id( e(2,id) ) = 0;
    end
end

Nf_TM = nnz(node_id);    % number of unknown nodes
Np_TM = Nn - Nf_TM;       % number of known nodes

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
    end
end


%{
    S_TE and T_TE are the stiffness and mass matrixes for TE modes,
    Sff_TM, Sfp_TM, Tff_TM, Tfp_TM, are the matrixes for TM modes.
    Ff_TM, Fp_TM, are the known and unknown nodes for TM modes.
%}
S_TE = spalloc(Nn,Nn,7*Nn);
T_TE = spalloc(Nn,Nn,7*Nn);
Sff_TM = spalloc(Nf_TM,Nf_TM,7*Nf_TM);
Tff_TM = spalloc(Nf_TM,Nf_TM,7*Nf_TM);
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
            Sij = (b(i)*b(j) + c(i)*c(j)) * Ae;

            % ---- TE modes ----
            S_TE(n(i),n(j)) = S_TE(n(i),n(j)) + Sij;

            if i == j
                T_TE(n(i),n(j)) = T_TE(n(i),n(j)) + Ae/6;
            else
                T_TE(n(i),n(j)) = T_TE(n(i),n(j)) + Ae/12;
            end

            % ---- TM modes ----
            if node_id(n(i)) == 1
                if node_id(n(j)) == 1
                    % both i and j are unknown
                    Sff_TM(index(n(i)),index(n(j))) = Sff_TM(index(n(i)),index(n(j))) + Sij;
                    if i == j
                        Tff_TM(index(n(i)),index(n(j))) = Tff_TM(index(n(i)),index(n(j))) + Ae/6;
                    else
                        Tff_TM(index(n(i)),index(n(j))) = Tff_TM(index(n(i)),index(n(j))) + Ae/12;
                    end
                end
            end
        end
    end
end

% Solve TE 
eig_num = 12;   % number of eigenvectors to get
[Hz_k_TE, ks_TE] = eigs(S_TE, T_TE, eig_num, "smallestabs");


% Solve TM
[Ef_z_k, ks_TM] = eigs(Sff_TM, Tff_TM, eig_num, "smallestabs");
% Now create the array for all nodes
Ez_k_TM = zeros(eig_num, Nn);
for mode = 1:eig_num
    for ind = 1:Nn
        if node_id(ind) == 1  
            % update only for unknown nodes
            Ez_k_TM(ind, mode) = Ef_z_k(index(ind), mode);
        end
    end
end


% Calculate the cutoff frequencies
fc_TE = zeros(1, eig_num);
for i = 1:eig_num
    fc_TE(i) = sqrt( ks_TE(i,i) ) * speed_of_light / (2 * pi);
end
fc_TM = zeros(1, eig_num);
for i = 1:eig_num
    fc_TM(i) = sqrt( ks_TM(i,i) ) * speed_of_light / (2 * pi);
end

disp(["The cutoff frequencies of the first ", eig_num, " calculated TE modes are:"])
disp(fc_TE);
disp(["The cutoff frequencies of the first ", eig_num, " calculated TM modes are:"])
disp(fc_TM);


% ------ calculate the x y fields ------

% TE
% the variables have an "a" infront of them to signify that they are 
% analogous to the actual values, but with different magnitude
aEx_k_TE = zeros(Ne, eig_num);
aEy_k_TE = zeros(Ne, eig_num);
aHx_k_TE = zeros(Ne, eig_num);
aHy_k_TE = zeros(Ne, eig_num);
for mode = 1:eig_num
    [ux, uy] = pdegrad(p,t,Hz_k_TE(:,mode));
    % Calculate the other coordinates without the constants 
    % because we don't have a frequency and eihter way we care
    % way we care about the flow 
    aEx_k_TE(:,mode) = - uy;
    aEy_k_TE(:,mode) =   ux;
    aHx_k_TE(:,mode) = - ux;
    aHy_k_TE(:,mode) = - uy;
end
% TM
aEx_k_TM = zeros(Ne, eig_num);
aEy_k_TM = zeros(Ne, eig_num);
aHx_k_TM = zeros(Ne, eig_num);
aHy_k_TM = zeros(Ne, eig_num);
for mode = 1:eig_num
    [ux, uy] = pdegrad(p,t,Ez_k_TM(:,mode));

    aEx_k_TM(:,mode) = - ux;
    aEy_k_TM(:,mode) = - uy;
    aHx_k_TM(:,mode) =   uy;
    aHy_k_TM(:,mode) = - ux;
end



% get the 9 modes with smallest cutoff frequency
first_9_modes = zeros(Nn, 9);

first_9_modes(:,1) = Ez_k_TM(:,3);
first_9_modes(:,2) = Hz_k_TE(:,1);
first_9_modes(:,3) = Ez_k_TM(:,4);
first_9_modes(:,4) = Hz_k_TE(:,3);
first_9_modes(:,5) = Ez_k_TM(:,6);
first_9_modes(:,6) = Ez_k_TM(:,8);
first_9_modes(:,7) = Hz_k_TE(:,4);
first_9_modes(:,8) = Ez_k_TM(:,10);
first_9_modes(:,9) = Ez_k_TM(:,12);








% ----------------- FIGURES ---------------------
warning('off', 'all');


% TE plot Hz 
figure('visible','off');
title("TE first " + eig_num + "modes");
for i = 1:eig_num
    subplot(ceil(sqrt(eig_num)), ceil(eig_num / ceil(sqrt(eig_num))), i);
    pdeplot(p, e, t, 'XYData', Hz_k_TE(:,i), Contour="on");
    title(['Mode ', num2str(i)]);
    colormap(jet);
    % Adjust colorbar precision
    cb = colorbar;
    ticks = cb.Ticks;
    cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'UniformOutput', false);
    axis equal tight;
end
exportgraphics(gcf, "./plots/TE_first_"+eig_num+"_modes.pdf", 'ContentType', 'vector');




% TM plot Ez
figure('visible','off');
title("TM first " + eig_num + "modes");
for i = 1:eig_num
    subplot(ceil(sqrt(eig_num)), ceil(eig_num / ceil(sqrt(eig_num))), i);
    pdeplot(p, e, t, 'XYData', Ez_k_TM(:,i), Contour="on");
    title(['Mode ', num2str(i)]);
    colormap(jet);
    axis equal tight;
end
exportgraphics(gcf, "./plots/TM_first_"+eig_num+"_modes.pdf", 'ContentType', 'vector');



set(gca, 'XTick', [], 'YTick', []);
% plot first nine modes 
figure('Visible','off');
layout = tiledlayout(3,3);
layout.TileSpacing = 'none';   
layout.Padding = 'none';  

% plot first nine modes 
figure('Visible','off');
layout = tiledlayout(3,3);
nexttile;
pdeplot(p,e,t, 'XYData', Hz_k_TE(:,3), Contour="on"); axis equal tight; colorbar;
title("TE"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Ez_k_TM(:,1), Contour="on"); axis equal tight; colorbar;
title("TM"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Hz_k_TE(:,4), Contour="on"); axis equal tight; colorbar;
title("TE"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Ez_k_TM(:,3), Contour="on"); axis equal tight; colorbar;
title("TM"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Hz_k_TE(:,6), Contour="on"); axis equal tight; colorbar;
title("TE"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Hz_k_TE(:,8), Contour="on"); axis equal tight; colorbar;
title("TE"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Ez_k_TM(:,4), Contour="on"); axis equal tight; colorbar;
title("TM"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Hz_k_TE(:,10), Contour="on"); axis equal tight; colorbar;
title("TE"); set(gca, 'XTick', [], 'YTick', []);
nexttile;
pdeplot(p,e,t, 'XYData', Hz_k_TE(:,12), Contour="on"); axis equal tight; colorbar;
title("TE"); set(gca, 'XTick', [], 'YTick', []);
colormap("jet");
exportgraphics(gcf, "./plots/Waveguide_first_nine_modes.pdf", 'ContentType', 'vector');




% plot the first 9 modes, transverse fields
figure('Visible','off');
layout = tiledlayout(3, 3);  
layout.TileSpacing = 'compact';   
layout.Padding = 'none';

nexttile;
plot_streamlines(aEx_k_TE(:,3), aEy_k_TE(:,3), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TE(:,3), aHy_k_TE(:,3), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TM(:,1), aEy_k_TM(:,1), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TM(:,1), aHy_k_TM(:,1), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TE(:,4), aEy_k_TE(:,4), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TE(:,4), aHy_k_TE(:,4), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TM(:,3), aEy_k_TM(:,3), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TM(:,3), aHy_k_TM(:,3), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TE(:,6), aEy_k_TE(:,6), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TE(:,6), aHy_k_TE(:,6), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TE(:,8), aEy_k_TE(:,8), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TE(:,8), aHy_k_TE(:,8), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TM(:,4), aEy_k_TM(:,4), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TM(:,4), aHy_k_TM(:,4), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TE(:,10), aEy_k_TE(:,10), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TE(:,10), aHy_k_TE(:,10), p, e, t, 'k'); axis equal tight off;

nexttile;
plot_streamlines(aEx_k_TE(:,12), aEy_k_TE(:,12), p, e, t, 'b'); hold on;
plot_streamlines(aHx_k_TE(:,12), aHy_k_TE(:,12), p, e, t, 'k'); axis equal tight off;

exportgraphics(gcf, "./plots/Waveguide_first_nine_modes_transverse.pdf", 'ContentType', 'vector');




%{
% TE plot transverse streamlines
figure('visible','off');
layout = tiledlayout(3, 3);
layout.TileSpacing = 'none';   
layout.Padding = 'none';       
for mode = 1:9
    [ux, uy] = pdegrad(p,t,Hz_k_TE(:,mode));
    aEx_TM = -uy;
    aEy_TM =  ux;
    aHx_TM = -ux;
    aHy_TM = -uy;

    nexttile;
    % Plot E field (blue)
    plot_streamlines(aEx_TM, aEy_TM, p, e, t, 'b');
    hold on;
    % Plot H field (red)
    plot_streamlines(aHx_TM, aHy_TM, p, e, t, 'k');
    axis tight;
    axis equal;
end
exportgraphics(gcf, "./plots/TE_first_"+eig_num+"_modes_transverse.pdf", 'ContentType', 'vector');

% TM plot transverse streamlines
figure('visible','off');
layout = tiledlayout(ceil(sqrt(eig_num)), ceil(eig_num / ceil(sqrt(eig_num))));
layout.TileSpacing = 'none';   
layout.Padding = 'none';     
for mode = 1:eig_num
    [ux, uy] = pdegrad(p,t,Ez_k_TM(:,mode));
    aEx_TM = -ux;
    aEy_TM = -uy;
    aHx_TM =  uy;
    aHy_TM = -ux;

    nexttile;
    % Plot E field (blue)
    plot_streamlines(aEx_TM, aEy_TM, p, e, t, 'b');
    hold on;
    % Plot H field (red)
    plot_streamlines(aHx_TM, aHy_TM, p, e, t, 'k');
    axis equal tight;
end
exportgraphics(gcf, "./plots/TM_first_"+eig_num+"_modes_transverse.pdf", 'ContentType', 'vector');
%}
