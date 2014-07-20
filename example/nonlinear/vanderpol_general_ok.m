% Vanderpol general 
function [ok] = vanderpol_general_ok(draw_plot)
    if nargin < 1
        draw_plot = false;
    end
    close all
    % ---------------------------------------------------------------
    % Simulation problem parameters
    max_itr = 30;           % Maximum number of optimizer iterations
    epsilon = 1e-4;         % Convergence criteria
    N       = 81;           % Number of measurement points (N > 1)
    T       = 8.;           % Total time
    ell     = 1;            % Number of constraints (never active)
    ndim    = 3;            % number of components tracked

    n       = 2*ndim;       % Number of components in the state vector (2 per component)
    m       = 1*ndim;       % Number of measurements at each time point 
 
    xi      = [ -.5 , -3., 0.3, .4, -0.2, .2]'; % True initial value for the state at time zero
    x_in    = zeros(n, N);  % Initial sequence where optimization begins
    initial = xi; % Estimate of initial state (at time index one)
    mu      = 2.0;          % ODE nonlinearity parameter
    sigma_r = .1;           % Standard deviation of measurement noise
    sigma_q = .1;           % Standard deviation of the transition noise
    sigma_e = .1;           % Standard deviation of initial state estimate
    randn('seed', 4321);    % Random number generator seed
    W      = 0.1*rand(2*ndim); % crosstalk matrix (random) 
    % -------------------------------------------------------------------
    % Level of tracing during optimization
    if draw_plot
        level = 1;
    else
        level = 1;
    end
    % spacing between time points
    dt = T / (N - 1);
    % -------------------------------------------------------------------
    %  information used by vanderpol_g

    params.a1 = 1; 
    params.a2 = 20; 
    params.a3 = 1; 
    params.a4 = 1; 
    params.a5 = 0;
    
    

    params.vanderpol_g_initial  = initial;
    params.vanderpol_g_dt       = dt;
    params.vanderpol_g_mu       = mu;
    params.W                    = W;
    % -------------------------------------------------------------------
    % global information used by direct_h
    %global direct_h_info
    params.direct_h_index = [1; 3; 5];
    
    % ---------------------------------------------------------------
    % Rest of the information required by ckbs_nonlinear
    %
    % Step size for fourth order Runge-Kutta method in vanderpol_sim
    % is the same as time between measurement points
    step    = dt;
    % Simulate the true values for the state vector
    x_true  = vanderpol_sim_full(params, xi, N, step, ndim, W);
    
       
    time    = linspace(0., T, N);
    % Simulate the measurement values
    H = [1 0 0 0 0 0; 0 0 1 0 0 0 ; 0 0 0 0 1 0]; 
    z       = H*x_true + sigma_r * randn(m, N);
    % Inverse covariance of the measurement and transition noise
    
    %
    rinv    = zeros(m, m, N);
    qinv    = zeros(n, n, N);
    for k = 1 : N
        rinv(:, :, k) = eye(m, m) / sigma_r^2;
        qinv(:, :, k) = eye(n, n) / sigma_q^2;
    end
    qinv(:, :, 1) = eye(n, n) / sigma_e^2;
    % ---------------------------------------------------------------
    % call back functions
    f_fun = @ no_f;
    g_fun = @ (k,x) vanderpol_full_g(k,x,params);
    h_fun = @(k,x) direct_h(k,x,params);
    % ---------------------------------------------------------------
    % call the optimizer
%
    [ x_out , u_out, info ] = ckbs_nonlinear( ...
        f_fun ,    ...
        g_fun ,    ...
        h_fun ,    ...
        max_itr ,  ...
        epsilon ,  ...
        x_in ,     ...
        z ,        ...
        qinv ,     ...
        rinv ,     ...
        level      ...
        );
    % ----------------------------------------------------------------------
    % Check that x_out is optimal
    % (this is an unconstrained case, so there is no need to check f).
%%
    ok     = size(info, 1) <= max_itr;
    g_out  = zeros(n,   N);
    h_out  = zeros(m,   N);
    dg_out = zeros(n, n,   N);
    dh_out = zeros(m, n,   N);
    xk     = zeros(n, 1);
    for k = 1 : N
        xkm      = xk;
        %
        xk       = x_out(:, k);
        uk       = u_out(:, k);
        [gk, Gk] = g_fun(k, xkm);
        [hk, Hk] = h_fun(k, xk);
        %
        g_out(:, k) = gk - xk;
        h_out(:, k) = hk;
        dg_out(:,:, k) = Gk;
        dh_out(:,:, k) = Hk;
    end
    dx    = zeros(n, N);
    grad  = ckbs_sumsq_grad(dx, z, g_out, h_out, dg_out, dh_out, qinv, rinv);
    ok    = max( max( abs(grad) ) ) < epsilon;
    % ----------------------------------------------------------------------
    if draw_plot
        figure(1);
        clf
        hold on
        plot(time, x_true(1,:), 'b-');
        plot(time, x_true(3,:)+1, 'b-');
        plot(time, x_true(5,:)-1, 'b-');

        plot(time, x_out(1,:),  'r--');
        plot(time, x_out(3,:)+1,  'r--');
        plot(time, x_out(5,:)-1,  'r--');

        plot(time, z(1,:),      'ko');
        plot(time, z(2,:)+1,      'k*');
        plot(time, z(3,:)-1,      'k+');

        title('Position: blue=truth, red=estimate');
        hold off
        
        figure(2)
        hold on
        plot(time, x_true(2,:), 'b-');
        plot(time, x_true(4,:)+2, 'b-');
        plot(time, x_true(6,:)-2, 'b-');

        plot(time, x_out(2,:),  'r--');
        plot(time, x_out(4,:)+2,  'r--');
        plot(time, x_out(6,:)-2,  'r--');
        title('Position: blue=truth, red=estimate');

        hold off
        savefig(1:2, ['/Users/saravkin/Dropbox/CompBioKalman/tex/figures/osc']);
       
        
        return
    end
end

 

% $$ $newlinech %$$
% $end
