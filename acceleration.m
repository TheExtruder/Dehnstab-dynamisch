function a = acceleration(M, K, u, F)

invM = inv(M);
a = -invM * (K * u - F);



