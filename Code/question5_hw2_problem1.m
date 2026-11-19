clear; clc; close all;

u_exact = @(x,y) cos(2*x).*sin(2*y);
f_fun   = @(x,y) -8*cos(2*x).*sin(2*y);
a = 0; b = 2*pi;

N = 50;
[x,y,U,Uexact,h] = solve_poisson(N,a,b,u_exact,f_fun);
[X,Y] = ndgrid(x,y);
% Colorblind-safe class palette (MATH550 textbook, "How to Visually
% Communicate"): unique, contrasting colors that stay distinguishable
% under protanopia and in grayscale.
cols = [ 67 179 174;   % teal          (marker edge/line color 1)
        136 216 192;   % light green   (marker face color 1)
        198  97 166;   % light magenta (marker face color 2)
        202  44 146]/255;   % magenta (marker edge/line color 2)

% 3D comparison of the exact and finite difference solutions.
fs3d = 35;  % 3D figure font size
figure('Color','w','Theme','light','Position',[100 100 1200 900]); hold on;
s1 = surf(X,Y,Uexact,'FaceColor',cols(1,:),'EdgeColor','none');
s2 = mesh(X,Y,U,'EdgeColor',cols(4,:),'FaceColor','none','LineWidth',0.8);
xlabel('$x$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
ylabel('$y$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
zlabel('$u(x,y)$','Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
title('HW2 Problem 1: Exact and Finite Difference Solutions', ...
    'Interpreter','latex','FontSize',0.8*fs3d,'Color','k');
legend([s1 s2],{'Exact solution','Finite difference solution'}, ...
    'Location','northeast','Interpreter','latex','FontSize',0.8*fs3d, ...
    'TextColor','black');
ax = gca;
ax.TickLabelInterpreter = 'latex';
% XColor/YColor/ZColor color the axis lines, tick marks, and tick labels.
set(ax,'LineWidth',3,'FontSize',0.8*fs3d, ...
    'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
ticks = [0 pi/2 pi 3*pi/2 2*pi];
tickLabel = {'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'};
ax.XTick = ticks;  ax.XTickLabel = tickLabel;
ax.YTick = ticks;  ax.YTickLabel = tickLabel;
view(45,30); grid on; box on; hold off;
exportgraphics(gcf,'poisson_hw2_problem1_solution.pdf','ContentType','vector');

Ns = [10 20 40 80 160];
errors = zeros(size(Ns)); hs = zeros(size(Ns));
for k = 1:numel(Ns)
    [~,~,Uh,Ue,hs(k)] = solve_poisson(Ns(k),a,b,u_exact,f_fun);
    errors(k) = max(abs(Uh(:)-Ue(:)));
end
orders = log(errors(1:end-1)./errors(2:end))./log(hs(1:end-1)./hs(2:end));
fprintf('    N           h           max error       observed order\n');
for k = 1:numel(Ns)
    if k == 1
        fprintf('%5d   %11.4e   %11.4e       ---\n',Ns(k),hs(k),errors(k));
    else
        fprintf('%5d   %11.4e   %11.4e       %.4f\n',Ns(k),hs(k),errors(k),orders(k-1));
    end
end

% Line-plot conventions (textbook): discrete data get symbols, theory/fits
% get curves, and reference lines are black, excluded from the legend, and
% labeled directly on the plot. Axes are log-log to reveal the power law.
fs = 50;  % convergence plot font size
figure('Color','w','Theme','light','Position',[100 100 1300 900]); hold on;
loglog(hs,errors(1)*(hs/hs(1)).^2,'--','Color','black','LineWidth',5, ...
    'HandleVisibility','off');  % O(h^2) reference: kept out of legend
loglog(hs,errors,'-d','Color',cols(1,:),'LineWidth',2,'MarkerSize',24, ...
    'MarkerEdgeColor',cols(1,:),'MarkerFaceColor',cols(2,:), ...
    'DisplayName','Numerical error');
hLab = 10^((log10(hs(1)) + log10(hs(end)))/2);  % geometric midpoint of h
text(1.5*hLab,0.4*errors(1)*(hLab/hs(1))^2,'$\mathcal{O}(h^2)$', ...
    'Interpreter','latex','FontSize',0.75*fs,'Color','k');
xlabel('$h$','Interpreter','latex','FontSize',fs,'Color','k');
ylabel('$\|u_h-u\|_\infty$','Interpreter','latex','FontSize',fs,'Color','k');
title('HW2 Problem 1: Poisson Convergence', ...
    'Interpreter','latex','FontSize',fs,'Color','k');
legend('Interpreter','latex','FontSize',0.75*fs,'Location','southeast', ...
    'TextColor','black');
ax = gca;
set(ax,'XScale','log','YScale','log', ...   % log-log error plot
    'LineWidth',3,'FontSize',fs,'TickLabelInterpreter','latex', ...
    'XColor',[0 0 0],'YColor',[0 0 0]);  % black axes, ticks, tick labels
grid on; box on; hold off;
exportgraphics(gcf,'poisson_hw2_problem1_convergence.pdf','ContentType','vector');

function [x,y,U,Uexact,h] = solve_poisson(N,a,b,u_exact,f_fun)
    h = (b-a)/N;
    x = linspace(a,b,N+1);
    y = linspace(a,b,N+1);
    xi = x(2:N); yi = y(2:N);
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
    U = u_exact_grid(x,y,u_exact);
    U(2:N,2:N) = reshape(L\F(:),n,n);
    Uexact = u_exact_grid(x,y,u_exact);
end

function U = u_exact_grid(x,y,u_exact)
    [X,Y] = ndgrid(x,y);
    U = u_exact(X,Y);
end
