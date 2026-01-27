
function x = PROXL21(U, lambda,tao)

lambda=tao*lambda;



U(U<0)=0;


row_norms = sqrt(sum(U.^2, 2));
scaling_factors = max(1 - lambda ./ row_norms, 0);
U = scaling_factors .* U;





x=U;
end

