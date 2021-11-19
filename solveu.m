function [u1,v1] = solveu(F, K, dirichletBoundary, M, u0, v0, dt)
%% Loesen des lineare Gleichungssystems unter Beruecksichtiung Dirichlet-Randbedingungen
%

invM = inv(M);
a = -invM * (K * u0 - F);
disp(a);
disp(u0);

u1 = zeros(size(K,1),1);         
u1 = u0 +dt * v0;
disp(u1);

disp(v0);    
v1 = zeros(size(K,1),1);         
v1= v0 + dt * a;
disp(v1);


end

