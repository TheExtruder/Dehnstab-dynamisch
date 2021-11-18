function Me = massmatrix(density, elementLength)

Me = zeros(2,2);
Me(1,1) = 2;     Me(1,2) = 1;

Me(2,1) = 1;     Me(2,2) = 2;

Me = Me * elementLength * density / 6;


