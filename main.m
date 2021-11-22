close all; clear all; clc;

disp('-------------------------------------------------------------------')
disp('-                           Dehnstab                              -') 
disp('-------------------------------------------------------------------')

%% Konstanten / Vorgaben

% Anzahl Elemente
disp('Geometry')
nElements = input('Enter number of elements: ');             % in x-Richtung
integrationParameter              = 1;           

% Geometrie
a   = 0;  b  = 4; 


    
% Materialdaten
emod    = 100;                 % E-Modul
density = 1;

% Materialmatrix
C = emod;

%% Netzgenerator
% Netz generieren (Elementfreiheitsgradzuordnungstabelle, Geometrievektor)

elementLength =  b / nElements;
[edof,q]       = mesh(a, nElements, elementLength);

% Geometriedatentabelle bzgl. Gebiet (1D Elemente)
[elementData]           = extract(edof,q);

% bzgl. Randpunkt
BoundaryEdofRight = edof(end , 2);
[BoundaryEdRight] = extract(BoundaryEdofRight,q);

%% Assembly
%   Allozieren / Speicher reservieren
nDof = size(q,1);
K    = sparse(nDof,nDof);
F    = sparse(nDof,1);
M    = sparse(nDof,nDof);

% Elementschleife
for i = 1:nElements
    [Ke,Fe] = element(elementLength,integrationParameter,C);
    [K]     = assem(K,Ke,edof(i,:));
    [F]     = assem(F,Fe,edof(i,:));
    
    Me = massmatrix(density, elementLength);
    [M]     = assem(M,Me,edof(i,:));
end

%% Schleife über Zeit




totalTime = 1;
dt  = 0.2;
u0  = zeros(size(K,1),1);
v0  = zeros(size(K,1),1);

for t = 0 : dt : totalTime

%   (horizontaler) Lastvektor
load = sin(t*pi); 

%% Randbedingungen einbauen
% Dirichlet-Rand
boundaryCondition1    = 0;
dirichletBoundary      = ones(1,2);
dirichletBoundary(1,2) = boundaryCondition1;
% Neumann-Rand
boundaryCondition2    = load;

% Randelement
F       = sparse(nDof,1);
[Fe]    = boundaryCondition2;
[F]     = assem(F,Fe(1),BoundaryEdofRight);



%% Massenmatrix


disp('t = ');
disp(t);
disp(u0);
[u1, v1] = solveu(F, K, M, u0, v0, dt, dirichletBoundary);
u0 = u1;
v0 = v1;
end
