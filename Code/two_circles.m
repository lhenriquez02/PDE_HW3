clear; clc; close all;

a = 0; b = 2*pi;
centers = [pi-1.25, pi; pi+1.25, pi];
radii = [0.65; 0.65];
gValues = [1; -1];
cols = [67 179 174; 136 216 192; 198 97 166; 202 44 146]/255;

N = 80;
[x,y,U,h,xb,yb,groups,gib,info] = solve_two_circles(N,a,b,centers,radii,gValues);
[X,Y] = ndgrid(x,y);
zib = gib;

fprintf('N=%d, h=%.4e, Nib=%d, GMRES flag=%d, relres=%.3e, true relres=%.3e, IB residual=%.3e\n', ...
    N,h,numel(xb),info.flag,info.relres,info.actual_relres,info.ib_residual);

% Divergent palette in the class colors, but with a dark center instead of
% white: teal -> charcoal -> magenta, so the u ~ 0 regions render dark
% (matching the example-code style) rather than washing out white.
% The endpoints are lifted toward white so they contrast against the
% dark charcoal center.
nHalf = 128;
tHalf = linspace(0,1,nHalf)';
darkMid = [24 26 32]/255;   % dark charcoal midpoint
lift   = 0.25;              % fraction of white blended into the endpoints
tealEnd = (1-lift)*cols(1,:) + lift*[1 1 1];
magEnd  = (1-lift)*cols(4,:) + lift*[1 1 1];
cmDiverging = [(1-tHalf)*tealEnd + tHalf*darkMid; ...
    (1-tHalf(2:end))*darkMid + tHalf(2:end)*magEnd];

% Symmetric color limits centered at zero.
cl = max(abs(U(:)));
if cl == 0, cl = 1; end

% Figure 1: 3D numerical solution and prescribed immersed boundary.
fs3d = 35;  % 3D solution figure font size
fig3d = figure('Color','w','Theme','light','Position',[100 100 1200 900]);
hold on;
s1 = surf(X,Y,U,'FaceAlpha',1,'EdgeColor','none');
colormap(cmDiverging); clim([-cl cl]);
p3 = gobjects(1,1);
for c = 1:numel(groups)
    ids = groups{c};
    hc = plot3([xb(ids);xb(ids(1))],[yb(ids);yb(ids(1))], ...
        [zib(ids);zib(ids(1))],'k-','LineWidth',2.5);
    if c == 1, p3 = hc; end
    plot3(xb(ids),yb(ids),zib(ids),'ko','MarkerSize',4, ...
        'MarkerFaceColor','k','HandleVisibility','off');
end
xlabel('$x$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
ylabel('$y$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
zlabel('$u(x,y)$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
title('HW2 Extension: Two Immersed Circles', ...
    'Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
legend([s1 p3],{'Finite difference solution','Immersed boundary'}, ...
    'Location','northeast','Interpreter','latex', ...
    'FontSize',0.7*fs3d,'TextColor','black');
ax = gca;
set(ax,'LineWidth',3,'FontSize',0.8*fs3d, ...
    'XColor','k','YColor','k','ZColor','k', ...
    'TickLabelInterpreter','latex');
ticks = [0 pi/2 pi 3*pi/2 2*pi];
tickLabel = {'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'};
ax.XTick = ticks; ax.XTickLabel = tickLabel;
ax.YTick = ticks; ax.YTickLabel = tickLabel;
cb = colorbar('TickLabelInterpreter','latex','FontSize',0.8*fs3d);
ylabel(cb,'$u(x,y)$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
view(45,30); grid on; box on; hold off;
exportgraphics(fig3d,'poisson_hw2_two_circles_solution.pdf','ContentType','vector');

% Figure 2: top view of the numerical solution, as in Problem 2.
fs2d = 35;  % top view figure font size
fig2d = figure('Color','w','Theme','light','Position',[100 100 1000 900]);
axT = axes; hold on;
daspect([1 1 1]);
pcolor(X,Y,U); shading interp;
colormap(fig2d,cmDiverging); clim([-cl cl]);
hPhi = plot(nan,nan,'s','MarkerSize',14, ...
    'MarkerFaceColor',tealEnd,'MarkerEdgeColor',tealEnd, ...
    'DisplayName','$u(x,y)$');
pIB = gobjects(1,1);
for c = 1:numel(groups)
    ids = groups{c};
    hc = plot([xb(ids);xb(ids(1))],[yb(ids);yb(ids(1))], ...
        'k-','LineWidth',2.5);
    if c == 1, pIB = hc; end
    plot(xb(ids),yb(ids),'ko','MarkerSize',4,'MarkerFaceColor','k', ...
        'HandleVisibility','off');
end
pIB.DisplayName = 'Immersed boundaries';
title('HW2 Extension: Two Immersed Circles (Top View)', ...
    'Interpreter','latex','FontSize',0.8*fs2d,'Color','k');
xlabel('$x$','Interpreter','latex','FontSize',0.8*fs2d,'Color','k');
ylabel('$y$','Interpreter','latex','FontSize',0.8*fs2d,'Color','k');
legend([hPhi pIB],'Location','northeast','Interpreter','latex', ...
    'FontSize',0.7*fs2d,'TextColor','black');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XTick = ticks; ax.XTickLabel = tickLabel;
ax.YTick = ticks; ax.YTickLabel = tickLabel;
set(ax,'LineWidth',3,'FontSize',0.8*fs2d,'XColor','k','YColor','k');
box on;
cb = colorbar(axT);
ylabel(cb,'$u(x,y)$','Interpreter','latex','FontSize',0.8*fs2d,'Color','k');
exportgraphics(fig2d,'poisson_hw2_two_circles_top.pdf','ContentType','vector');

function [x,y,U,h,xb,yb,groups,g,info] = ...
    solve_two_circles(N,a,b,centers,radii,gValues)

    h = (b-a)/N;
    x = linspace(a,b,N+1);
    y = linspace(a,b,N+1);

    xi = x(2:N);
    yi = y(2:N);

    n = N-1;
    e = ones(n,1);

    D = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;
    I = speye(n);
    L = kron(I,D) + kron(D,I);

    F = zeros(n,n);

    xb = []; yb = []; g = [];
    groups = cell(numel(radii),1);
    for c = 1:numel(radii)
        nMarkers = max(8,ceil(2*pi*radii(c)/h));
        theta = (0:nMarkers-1)'*(2*pi/nMarkers);
        ids = numel(xb) + (1:nMarkers);
        groups{c} = ids(:);
        xb = [xb; centers(c,1)+radii(c)*cos(theta)];
        yb = [yb; centers(c,2)+radii(c)*sin(theta)];
        g = [g; gValues(c)*ones(nMarkers,1)];
    end
    Nib = numel(xb);

    gaussian_width = 1.2*h;

    [S,J,~] = ib_operators(xi,yi,xb,yb,h,gaussian_width);

    Uboundary = zeros(N+1,N+1);
    gb = g;

    rhs = [F(:); gb];
    nu = n*n;
    nq = Nib;

    Amult = @(z) [L*z(1:nu) - S*z(nu+1:end); ...
                   J*z(1:nu)];

    Lsolve = decomposition(L,'lu');
    T = Lsolve\S;
    Schur = full(J*T);
    SchurSolve = decomposition(Schur,'lu');

    Psolve = @(r) block_precondition(r,Lsolve,T,SchurSolve,nu);

    tol = 1e-9;
    restart = min(40,nu+nq);
    maxit = 5;
    [sol,flag,relres,iter] = gmres( ...
        Amult,rhs,restart,tol,maxit,Psolve);

    actual_relres = norm(Amult(sol)-rhs)/max(norm(rhs),eps);
    if flag ~= 0 || actual_relres > 1e-7
        warning('N=%d: GMRES flag=%d, relres=%.3e, true relres=%.3e.', ...
            N,flag,relres,actual_relres);
    end

    U = Uboundary;
    U(2:N,2:N) = reshape(sol(1:nu),n,n);

    info.flag = flag;
    info.relres = relres;
    info.iterations = (iter(1)-1)*restart+iter(2);
    info.ib_residual = norm(J*sol(1:nu)-gb,inf);
    info.actual_relres = actual_relres;
end

function [S,J,B] = ib_operators(xi,yi,xb,yb,h,a)

    n = numel(xi);
    N = n+1;
    fullx = [xi(1)-h, xi(:).', xi(end)+h];
    fully = [yi(1)-h, yi(:).', yi(end)+h];
    nq = numel(xb);
    radius = 8*a;

    interior_rows = cell(nq,1);
    interior_cols = cell(nq,1);
    interior_vals = cell(nq,1);
    boundary_rows = cell(nq,1);
    boundary_cols = cell(nq,1);
    boundary_vals = cell(nq,1);

    for k = 1:nq
        ix = find(abs(fullx-xb(k)) <= radius);
        iy = find(abs(fully-yb(k)) <= radius);
        [II,JJ] = ndgrid(ix,iy);
        II = II(:); JJ = JJ(:);
        r2 = (fullx(II)-xb(k)).^2 + (fully(JJ)-yb(k)).^2;
        r2 = r2(:);
        keep = r2 <= radius^2;
        II = II(keep); JJ = JJ(keep); r2 = r2(keep);
        weights = exp(-r2/(2*a^2))/(2*pi*a^2);

        isInterior = II>=2 & II<=N & JJ>=2 & JJ<=N;
        Iin = II(isInterior); Jin = JJ(isInterior);
        interior_rows{k} = (Jin-2)*n + (Iin-1);
        interior_cols{k} = repmat(k,numel(Iin),1);
        interior_vals{k} = weights(isInterior);

        Ibd = II(~isInterior); Jbd = JJ(~isInterior);
        boundary_rows{k} = repmat(k,numel(Ibd),1);
        boundary_cols{k} = (Jbd-1)*(N+1) + Ibd;
        boundary_vals{k} = h^2*weights(~isInterior);
    end

    W = sparse(vertcat(interior_rows{:}),vertcat(interior_cols{:}), ...
               vertcat(interior_vals{:}),n*n,nq);
    S = W;
    J = h^2*W';
    B = sparse(vertcat(boundary_rows{:}),vertcat(boundary_cols{:}), ...
               vertcat(boundary_vals{:}),nq,(N+1)^2);
end

function z = block_precondition(r,Lsolve,T,SchurSolve,nu)
    zq = SchurSolve\r(nu+1:end);
    zu = Lsolve\r(1:nu) + T*zq;
    z = [zu;zq];
end
