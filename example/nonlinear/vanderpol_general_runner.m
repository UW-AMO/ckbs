k = 2; % say we are at a generic point; 
ndim = 3; % number of components (modeled via vanderpol)
W = rand(ndim*2); % cross talk matrix linking components

initial = rand(ndim*2, 1); 
vanderpol_g_dt = 0.01;
vanderpol_g_mu = 2.0; 

params.vanderpol_g_initial = initial; % should be a vector of size 2n.
params.vanderpol_g_dt = vanderpol_g_dt;
params.vanderpol_g_mu = vanderpol_g_mu;

params.W = W; 
params.ndim = ndim;


xk1 = rand(ndim*2, 1); 

[gk, Gk] = vanderpol_general_g(k, xk1, params);


