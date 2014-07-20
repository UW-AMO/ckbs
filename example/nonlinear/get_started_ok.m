% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin get_started_ok.m$$ $newlinech %$$
% $spell
%	matlab
%	ko
%	complementarity
%	nargin
%	randn
%	optimizer
%	clf
%	gk
%	est
%	mu
%	itr
%	uk
%	dk
%	ckbs
%	qinv
%	qinvk
%	rinv
%	rinvk
%	xk
%	ek
%	df
%	dg
%	dh
%	Fk
%	Gk
%	Hk
%	sumsq
%   params
% $$
%
% $section ckbs_nonlinear: A Simple Example and Test$$
%
% $index ckbs_nonlinear, simple$$
% $index ckbs_nonlinear, get_started$$
% $index example, get_started$$
% $index test, get_started$$
%
% $index get_started, ckbs_nonlinear$$
% $index simple, ckbs_nonlinear$$
% $index example, ckbs_nonlinear$$
% $index test, ckbs_nonlinear$$
%
% $head Syntax$$
% $codei%[%ok%] = get_started_ok(%draw_plot%)%$$
%
% $head Running This Example$$
% Change into the $code example$$ directory,
% start $code octave$$ (or $code matlab$$), and then
% execute the commands
% $codep
%	test_path
%	get_started_ok(true)
% %$$
%
% $head draw_plot$$
% If this optional argument is true, a plot is drawn displaying the results.
% If $icode draw_plot$$ is no present or false, no plot is drawn. 
%
% $head ok$$
% If the return value $icode ok$$ is true, the test passed,
% otherwise the test failed.
%
% $head Simulated State Vector$$
% The state vector for this example is in $latex \B{R}^n$$ 
% and the sequence of state vector values is simulated with 
% $latex x_{j,k} = 1$$ 
% for $latex j = 1 , \ldots , n$$ and $latex k = 1 , \ldots N$$.
%
% $head Transition Model$$
% We simulate imperfect knowledge of this dynamical system
% using the persistence model. 
% To be specific,
% $latex \[
%	x_k = x_{k-1} + w_k
% \] $$ 
% where $latex w_k \sim \B{N} ( 0 , Q_k )$$ 
% and the variance $latex Q_k = I_n$$ (and $latex I_n$$ is the
% $latex n \times n$$ identity matrix.
% The routine $cref persist_g.m$$ calculates the mean for $latex x_k$$
% given $latex x_{k-1}$$.
% Note that the initial state estimate is returned by this routine as
% $latex \hat{x} = g_1 (x)$$ and has $latex \hat{x}_{j,k} = 1$$
% for $latex j = 1 , \ldots , n$$.
% Also not that the corresponding variance $latex Q_1 = I_n$$. 
%
% $head Measurements$$
% For this example, the measurements are direct observations of the state
% $latex \[
%	z_k = x_k + v_k
% \] $$
% where $latex v_k \sim \B{N} ( 0 , R_k )$$ 
% and the variance $latex R_k = I_n$$.
% The routine $cref direct_h.m$$ calculates the mean for $latex z_k$$,
% given the value of $latex x_k$$.
%
% $head Constraints$$
% For this example, the constraints are 
% $latex \[
%	0.5 \leq x_{k,j} \leq 1.5
% \] $$
% for $latex j = 1 , \ldots , n$$ and $latex k = 1 , \ldots , N$$.
% The routine $cref box_f.m$$ represents these constraints in the form
% $latex f_k ( x_k ) \leq 0$$ as expected by $cref ckbs_nonlinear$$.
%	
% $head Call Back Functions$$
% $table
% $rref persist_g.m$$
% $rref direct_h.m$$
% $rref box_f.m$$
% $tend
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = get_started_ok(draw_plot)
    if nargin < 1
        draw_plot = false;
    end
    % --------------------------------------------------------------------------
    max_itr = 20;      % maximum number of iterations
    epsilon = 1e-5;    % convergence criteria
    N       = 40;      % number of time points
    n       = 1;       % number of state vector components per time point
    sigma_r = 1;       % variance of measurement noise
    sigma_q = 1;       % variance of transition noise
    sigma_e = 1;       % variance of initial state estimate
    initial = ones(n); % estimate for initial state value
    randn('seed', 4321)
    % ---------------------------------------------------------------------------
    % Level of tracing during optimization
    if draw_plot
        level = 1;
    else
        level = 0;
    end
    m       = n;       % number of measurements per time point
    ell     = 2 * n;   % number of constraints
    index   = 1 : n;   % index vector
    % ---------------------------------------------------------------------------
    %  information used by box_f
    params.box_f_lower = 0.5 * ones(n, 1);
    params.box_f_upper = 1.5 * ones(n, 1);
    params.box_f_index = index;
    % ---------------------------------------------------------------------------
    % global information used by persist_g
    params.persist_g_initial = initial;
    % ---------------------------------------------------------------------------
    % global information used by direct_h
    %global direct_h_info
    
    params.direct_h_index = index;
     h_fun = @(k,x) direct_h(k,x,params);
    % ---------------------------------------------------------------------------
    % simulated true trajectory 
    x_true  = ones(n, N);
    % initial vector during optimization process
    x_in    = 0.5 * ones(n, N);
    % measurement variances
    rinv    = zeros(m, m, N);
    qinv    = zeros(n, n, N);
    z       = zeros(m, N);
    for k = 1 : N
        xk           = x_true(:, k);
        ek           = randn(m, 1);
        hk           = h_fun(k, xk);
        z(:, k)      = hk + ek;
        qinv(:,:, k) = eye(n) / (sigma_q * sigma_q);
        rinv(:,:, k) = eye(n) / (sigma_r * sigma_r);
    end
    qinv(:,:,1)       = eye(n) / (sigma_e * sigma_e);
    % ----------------------------------------------------------------------
    % call back functions
    f_fun = @(k,xk) box_f(k,xk,params);
    g_fun = @(k,xk) persist_g(k,xk,params);
   
    % ----------------------------------------------------------------------
    % call optimizer
    [x_out, u_out, info] = ckbs_nonlinear( ...
        f_fun,    ...
        g_fun,    ...
        h_fun,    ...
        max_itr,  ...
        epsilon,  ...
        x_in,     ...
        z,        ...
        qinv,     ...
        rinv,     ...
        level     ...
    );
    % ----------------------------------------------------------------------
    ok     = size(info, 1) <= max_itr;
    %
    % Set up affine problem corresponding to x_out (for a change in x).
    % Check constraint feasibility and complementarity.
    f_out  = zeros(ell, N);
    g_out  = zeros(n,   N);
    h_out  = zeros(m,   N);
    df_out = zeros(ell, n, N);
    dg_out = zeros(n, n,   N);
    dh_out = zeros(m, n,   N);
    xk1    = zeros(n, 1);
    for k = 1 : N
        xk       = x_out(:, k);
        uk       = u_out(:, k);
        [fk, Fk] = f_fun(k, xk);
        [gk, Gk] = g_fun(k, xk1);
        [hk, Hk] = h_fun(k, xk);
        %
        ok       = ok & all( fk <= epsilon);              % feasibility
        ok       = ok & all( abs(fk) .* uk <= epsilon );  % complementarity
        %
        f_out(:, k)    = fk;
        g_out(:, k)    = gk - xk;
        h_out(:, k)    = hk;
        df_out(:,:, k) = Fk;
        dg_out(:,:, k) = Gk;
        dh_out(:,:, k) = Hk;
        xk1 = xk;
    end
    ok     = ok & all( all( u_out >= 0 ) );
    %
    % Compute gradient of unconstrained affine problem at dx equal to zero
    % and then check the gradient of the Lagrangian
    dx_out = zeros(n, N);
    grad  = ckbs_sumsq_grad(dx_out, z, g_out, h_out, dg_out, dh_out, qinv, rinv);
    for k = 1 : N
        uk  = u_out(:, k);
        Fk  = df_out(:,:, k);
        dk  = grad(:, k);
        %
        ok  = ok & all( abs(Fk' * uk + dk) <= epsilon );
    end
    % ----------------------------------------------------------------------
    if draw_plot
        figure(1);
        clf
        hold on
        time = 1 : N;
        plot(time, x_true(1,:),    'b-');
        plot(time, x_out(1,:),     'g-');
        plot(time, z(1,:),         'ko'); 
        plot(time, 0.5*ones(1,N),  'r-');
        plot(time, 1.5*ones(1,N),  'r-');
        title('Position: blue=truth, green=estimate, red=constraint, o=data');
        axis( [ 0 , N , -1, 3 ] )
        hold off
    return
end
% $$ $newlinech %$$
% $end
