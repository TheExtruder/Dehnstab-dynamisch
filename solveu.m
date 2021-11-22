function [u1,v1] = solveu(F, K, M, u0, v0, dt, dirichletBoundary)
%% Loesen des lineare Gleichungssystems unter Beruecksichtiung Dirichlet-Randbedingungen

if ~isempty(dirichletBoundary)
    % Eingaben einlesen
    dofDiri = dirichletBoundary(:,1);   % 1
    qDiri   = dirichletBoundary(:,2);   % 0

    % Freiheitsgrade mit Dirichlet-Raender identifizieren
    nDof = numel(F);
    DoSolve = true(nDof,1);         % logisches Array
    DoSolve(dofDiri) = false;       % erster Eintrag = 0

    u1 = zeros(size(K,1),1); 
    v1 = zeros(size(K,1),1); 
    invM = inv(M);

    a = -invM * (K * u0 - F);
    
    u1(DoSolve) = u0(DoSolve) + dt * v0(DoSolve);        
    u1(~DoSolve) = qDiri;   

    v1 = v0 + dt * a;

else
disp('Error in solve function')
end
end

