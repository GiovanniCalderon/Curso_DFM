
% Programa principal para resolver la ecuación: 
% 
%    -(k u_x)_x = 0 en [0,1],  k(x)=1 en el intervalo (0,1)
% Con condiciones de Dirchlet :
%       u(0) = 0,  u(1) = 1
%  
% con  un esquema mimético.
% PARAMETROS NECESARIOS :
%    nelem : Número de elementos de la malla
%    [a,b] : Valores extremos del dominio de definicion
%    fuente: Nombre de la funcion que define el termino fuente del problema
% Observación: 
%   . Se usa la libreria MOLE 
%   . La solución exacta esta dada por: 
%         u(x)=x,  q(x)=-1; 
%  El código esta pensando para introducir los Miméticos en un curso corto
%  usando la libreria MOLE  
%                GECS   11/10/2025
%-------------------------------------------------------------

%% Ejemplo MDM/MOLE
%  -(k u_x)_x = 0 en [0,1],  k(x)=1
%  u(0) = 0,  u(1) = 1
%  Solución exacta: u(x)=x,  q(x)=-1

clear; clc;
%% Añadir MOLE al path (ajusta la ruta en tu máquina si hace falta)
addpath('C:\Users\TATUY\OneDrive - Universidad Industrial de Santander\INVESTIGACION_UIS\Curso_MDM_MOLE\mole-main\mole-main\src\matlab_octave')

%% Parámetros del problema
west  = 0;
east  = 1;
orden = 2;     % orden mimético
nelem = 50;    % número de celdas
nnode = nelem+1;  % number de nodos (enteros) para MOLE
nnode2= nnode+1;  % number de nodos (partidos)+frontera
dx = (east - west)/nelem; % tamaño de la celda

%% Operadores miméticos 1D
G = grad(orden, nelem, dx);    % gradiente (n_faces x n_nodes)
D = div(orden, nelem, dx);     % divergencia (n_nodes x n_faces)

%% Geometría (uniforme)
x_nodes = linspace(west,east,nnode)';             % nodos enteros
x_centro = [west west+.5*dx: dx :east-.5*dx east]; % puntos medios

%% Conductividad en las caras (constante k=1)
k_nodes = ones(nnode,1);
Kf = spdiags(k_nodes(:),0,nnode,nnode);        % matriz diagonal

%% Operador elíptico mimético
A = -D * Kf * G;           % A u ≈ -(k u_x)_x

%% Término fuente (f(x)=0)
b = zeros(nnode2,1);

%% Condiciones de frontera (Robin general a*u + b*u_x = g)
% Aquí: Dirichlet-Dirichlet
% u(0)=0  -> a0=1, b0=0, g0=0
% u(1)=1  -> aL=1, bL=0, gL=1

dc = [1; 1];      % a en x=0 y x=L
nc = [0; 0];      % b en x=0 y x=L
vc = [0; 1];      % g en x=0 y x=L

% Imponer BC con la rutina de MOLE
[A_bc, b_bc] = addScalarBC1D(A, b, orden, nelem, dx, dc, nc, vc);

%% Resolver sistema lineal
u_h = A_bc \ b_bc;       % solución mimética en los nodos

%% Solución exacta y error
u_exact = x_nodes;       % u(x)=x
err_inf = norm(u_h' - u_exact, inf);

fprintf('Error infinito ||u_h - u||_inf = %.3e\n', err_inf);

%% Flujo numérico q_h = -k u_x en las caras
grad_u = G * u_h;              % u_x en caras
q_h    = -k_nodes .* grad_u;   % q = -k u_x

fprintf('q_h(1) = %.6f, q_h(end) = %.6f (exacto: -1)\n', ...
         q_h(1), q_h(end));

%% Gráficas
figure;
plot(x_centro, u_h, '-o', x_nodes, u_exact, '--');
xlabel('x'); ylabel('u(x)');
legend('u_h (MDM)','u(x)=x','Location','Best');
title('Solución MDM/MOLE para -(k u_x)_x = 0, k=1');
grid on;

figure;
plot(x_nodes, q_h, '-s');
xlabel('x'); ylabel('q(x)');
title('Flujo mimético q_h(x) = -k u_x');
grid on;
