
clear; clc; close all;

u_exact = @(x,y) cos(2*x).*sin(2*y);
f_fun   = @(x,y) -8*cos(2*x).*sin(2*y);
a = 0; b = 2*pi;

cols = [67 179 174; 136 216 192; 198 97 166; 202 44 146]/255;

N = 50;
[x,y,U,Uexact,h,xb,yb,info] = solve_poisson_ib(N,a,b,u_exact,f_fun);
[X,Y] = ndgrid(x,y);

fprintf('N=%d, h=%.4e, Nib=%d, GMRES flag=%d, relres=%.3e, iterations=%d\n', ...
    N,h,numel(xb),info.flag,info.relres,info.iterations);

% 3D comparison of the exact and immersed boundary solutions.
% Field-plot convention: the oscillatory field gets a divergent colormap
% built from the class palette (teal -> white -> magenta), so different
% heights carry different colors, centered on zero.
fs3d = 35;  % 3D solution figure font size

nHalf = 128;
tHalf = linspace(0,1,nHalf)';
cmDiverging = [(1 - tHalf)*cols(1,:) + tHalf*[1 1 1]; ...
    (1 - tHalf(2:end))*[1 1 1] + tHalf(2:end)*cols(4,:)];

figure('Color','w','Theme','light','Position',[100 100 1200 900]); hold on;

s1 = surf(X,Y,Uexact,'EdgeColor','none');
colormap(cmDiverging);
clim([-1 1]);     % center the divergent map on zero height

s2 = mesh(X,Y,U,'EdgeColor',0.6*cols(4,:), ...   % darker magenta: stands off the magenta peaks
    'FaceColor','none','LineWidth',0.8);

% Plot the immersed curve at the prescribed boundary height.
zib = u_exact(xb,yb);
p3 = plot3([xb;xb(1)],[yb;yb(1)],[zib;zib(1)], ...
    'k-','LineWidth',2.5);

plot3(xb,yb,zib,'ko','MarkerSize',4, ...
    'MarkerFaceColor','k','HandleVisibility','off');

xlabel('$x$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
ylabel('$y$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
zlabel('$u(x,y)$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');

title('HW2 Problem 2: Exact and Immersed Boundary Solutions', ...
    'Interpreter','latex','FontSize',0.8*fs3d,'Color','k');

legend([s1 s2 p3], ...
    {'Exact solution','Finite difference solution','Immersed boundary'}, ...
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

exportgraphics(gcf,'poisson_hw2_problem2_solution.pdf', ...
    'ContentType','vector');

% (handout Fig. 2 / example-code style): solution with a divergent
% colormap from the class palette, black IB curve with markers, legend.
fs2d = 35;
fig2d = figure('Color','w','Theme','light','Position',[100 100 1000 900]);
axT = axes; hold on;
daspect([1 1 1]);

% Approximate solution as the flat colored field.
pcolor(X,Y,U);
shading interp;
colormap(fig2d,cmDiverging);
clim([-1 1]);

% Legend swatch for the field (flat color square, invisible on plot).
hPhi = plot(nan,nan,'s','MarkerSize',14, ...
    'MarkerFaceColor',cols(1,:),'MarkerEdgeColor',cols(1,:), ...
    'DisplayName','$u(x,y)$');

% Immersed boundary curve and markers.
pIB = plot([xb;xb(1)],[yb;yb(1)],'k-','LineWidth',2.5);
pIB.DisplayName = 'Immersed boundary';
plot(xb,yb,'ko','MarkerSize',4,'MarkerFaceColor','k', ...
    'HandleVisibility','off');

title('HW2 Problem 2: Finite Difference Solution (Top View)', ...
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

exportgraphics(fig2d,'poisson_hw2_problem2_top.pdf', ...
    'ContentType','vector');

Ns = [10 20 40 80 160 320];
errors = zeros(size(Ns));
errL2   = zeros(size(Ns));
hs = zeros(size(Ns));

for k = 1:numel(Ns)
    [~,~,Uh,Ue,hs(k),~,~,stats] = ...
        solve_poisson_ib(Ns(k),a,b,u_exact,f_fun);

    errors(k) = max(abs(Uh(:)-Ue(:)));
    errL2(k)  = norm(Uh(:)-Ue(:))*hs(k);

    fprintf('N=%d: GMRES flag=%d, relres=%.3e\n', ...
        Ns(k),stats.flag,stats.relres);
end

orders = log(errors(1:end-1)./errors(2:end))./ ...
         log(hs(1:end-1)./hs(2:end));
ordersL2 = log(errL2(1:end-1)./errL2(2:end))./ ...
           log(hs(1:end-1)./hs(2:end));

fprintf('\n    N        L-inf error    order       L2 error     order\n');
for k = 1:numel(Ns)
    if k == 1
        fprintf('%5d   %11.4e   ------    %11.4e   ------\n', ...
            Ns(k),errors(k),errL2(k));
    else
        fprintf('%5d   %11.4e   %6.3f    %11.4e   %6.3f\n', ...
            Ns(k),errors(k),orders(k-1),errL2(k),ordersL2(k-1));
    end
end

fs = 40;  % Convergence plot font size
figure('Color','w','Theme','light','Position',[100 100 1300 900]); hold on;

% Reference slope anchored at the finest-grid error: the asymptotic
% segments then overlay the line, making the h^1.5 rate easy to read.
loglog(hs,errL2(end)*(hs/hs(end)).^1.5,'--', ...
    'Color','k','LineWidth',4,'HandleVisibility','off');

loglog(hs,errL2,'-d','Color',cols(1,:), ...
    'LineWidth',2,'MarkerSize',20, ...
    'MarkerEdgeColor',cols(1,:), ...
    'MarkerFaceColor',cols(2,:), ...
    'DisplayName','Numerical error');

hLab = sqrt(hs(1)*hs(end));

text(1.25*hLab,0.5*errL2(end)*(hLab/hs(end))^1.5, ...
    '$\mathcal{O}(h^{3/2})$', ...
    'Interpreter','latex','FontSize',0.75*fs,'Color','k');

xlabel('$h$','Interpreter','latex','FontSize',fs,'Color','k');
ylabel('$\|u_h-u\|_{2}$', ...
    'Interpreter','latex','FontSize',fs,'Color','k');

title('HW2 Problem 2: Immersed Boundary Convergence', ...
    'Interpreter','latex','FontSize',fs,'Color','k');

legend('Interpreter','latex','FontSize',0.75*fs, ...
    'Location','southeast','TextColor','black');

ax = gca;
set(ax,'XScale','log','YScale','log', ...
    'LineWidth',3,'FontSize',fs, ...
    'TickLabelInterpreter','latex', ...
    'XColor','k','YColor','k');

grid on; box on; hold off;

exportgraphics(gcf,'poisson_hw2_problem2_convergence.pdf', ...
    'ContentType','vector');

function [x,y,U,Uexact,h,xb,yb,info] = ...
    solve_poisson_ib(N,a,b,u_exact,f_fun)

    h = (b-a)/N;
    x = linspace(a,b,N+1);
    y = linspace(a,b,N+1);

    xi = x(2:N);
    yi = y(2:N);

    [X,Y] = ndgrid(xi,yi);

    n = N-1;
    e = ones(n,1);

    D = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;
    I = speye(n);
    L = kron(I,D) + kron(D,I);

    F = f_fun(X,Y);

    F(1,:)   = F(1,:)   - u_exact(a,yi)/h^2;
    F(end,:) = F(end,:) - u_exact(b,yi)/h^2;
    F(:,1)   = F(:,1)   - u_exact(xi,a)'/h^2;
    F(:,end) = F(:,end) - u_exact(xi,b)'/h^2;

    Nib = max(8,round(2*pi/h));

    theta = (0:Nib-1)'*(2*pi/Nib);

    xb = pi + cos(theta);
    yb = pi + sin(theta);

    ds = 2*pi/Nib;
    g = u_exact(xb,yb);

    gaussian_width = 1.2*h;

    [S,J,B] = ib_operators(xi,yi,xb,yb,h,gaussian_width);

    Uboundary = u_exact_grid(x,y,u_exact);
    Uboundary(2:N,2:N) = 0;
    gb = g - B*Uboundary(:);

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

    Uexact = u_exact_grid(x,y,u_exact);

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

function U = u_exact_grid(x,y,u_exact)
    [X,Y] = ndgrid(x,y);
    U = u_exact(X,Y);
end
