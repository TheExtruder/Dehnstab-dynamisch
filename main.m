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

%   (horizontaler) Lastvektor
load = 25; 
    
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
end

%% Randbedingungen einbauen
% Dirichlet-Rand
boundaryCondition1    = 0;
dirichletBoundary      = ones(1,2);
dirichletBoundary(1,2) = boundaryCondition1;
% Neumann-Rand
boundaryCondition2    = load;

% Randelement
    [Fe]    = boundaryCondition2;
    [F]     = assem(F,Fe(1),BoundaryEdofRight);


%% Gleichungsloeser
d = solveq(K,F,dirichletBoundary);

%% Darstellung Ausgangssituation
y = zeros(1, nElements + 1);
plot(q,y,'r-o','LineWidth',3);
hold on

%% Darstellung Gleichgewichtslage
q = q + d;                      % Update, Konfiguration Gleichgewichtslage
plot(q,y,'b-o','LineWidth',2);
hold off

[newElementData] = extract(edof,q); 

disp('Ausgangslage:'  )
disp(elementData);
disp('Verschiebungsvektor:'  )
disp(d);
disp('Gleichgewichtslage:'  )
disp(newElementData);

%% Massenmatrix

Me = massmatrix(density, elementLength);
for i = 1:nElements
    [M]     = assem(M,Me,edof(i,:));
end
invM = inv(M);
a = -invM * (K * d - load); %acceleration
disp(a);

dt  = 1;
dn  = d;
vn  = 1;

dn1 = dn + dt * vn;
vn1 = vn + dt * a;

dn2 = dn1 + dt * vn;
vn2 = vn1 + dt * a;

dt = 1;

