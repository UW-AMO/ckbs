% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2010
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          James Burks:          burke at math dot washington dot edu
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin vanderpol_ok.m$$ $newlinech %$$
% $spell
%       ko
%       hk
%       Hk
%       sumsq
%       clf
%       src
%       dg
%       dh
%       xkm
%       uk
%       gk
%       Gk
%       linspace
%       randn
%       rinv
%       qinv
%       nonlinearity
%       Runge-Kutta
%       sim
%       ckbs
%       mu
%       est
%       itr
%       xi
%       optimizer
%       dt
%       xk
%       tk
%       dxk
%       nargin
%       vanderpol
%       Gillijns
%       Kandepu Foss Imsland
%       params
%       end end
% $$
%
% $section ckbs_nonlinear:
% Unconstrained Nonlinear Transition Model Example and Test$$
%
% $index ckbs_nonlinear, unconstrained$$
% $index example, unconstrained$$
% $index test, unconstrained$$
%
% $index vanderpol, unconstrained$$
% $index unconstrained, ckbs_nonlinear$$
%
% $index example, unconstrained$$
% $index test, unconstrained$$
%
% $index example, transition nonlinear$$
% $index test, transition nonlinear$$
%
% $head Syntax$$
% $codei%[%ok%] = vanderpol_ok(%draw_plot%)%$$
%
% $head draw_plot$$
% If this optional argument is true, a plot is drawn showing the results.
% If $icode draw_plot$$ is not present or false, no plot is drawn.
%
% $head ok$$
% If the return value $icode ok$$ is true, the test passed,
% otherwise the test failed.
%
% $head Reference$$
% $index reference$$
% This example appears in many places; for example
%
% $list number$$
% Section 4.1.1 of
% $italic
% Applying the unscented Kalman filter for nonlinear state estimation
% $$, R. Kandepu, B. Foss , L. Imsland,
% Journal of Process Control,
% vol. 18, pp. 753, 2008.
%
% $lnext
% Section V of
% $italic
% What Is the Ensemble Kalman Filter and How Well Does it Work?
% $$,
% S. Gillijns et-al.
% Proceedings of the 2006 American Control Conference,
% Minneapolis, June 14-16, 2006
%
% $lend
%
% $head Simulated State Vector$$
% We use $latex x1 (t)$$ and $latex x2 (t)$$
% to denote the oscillator position and velocity as a function of time.
% The deterministic ordinary differential equation for the
% Van der Pol oscillator is
% $latex \[
% \begin{array}{rcl}
%       x1 '(t) & = & x2 (t)
%       \\
%       x2 '(t) & = & \mu [ 1 - x1(t)^2 ] x2 (t) - x1(t)
% \end{array}
% \] $$
% The simulated state vector values are the solution to this ODE
% with nonlinear parameter $latex \mu = 2$$ and initial conditions
% $latex x1 (0) = 0$$, $latex x2 (0) = -5.$$.
% The routine $cref vanderpol_sim$$ calculates these values.
%
% $head Transition Model$$
% We simulate imperfect knowledge in the dynamical system
% by using Euler's approximation for the ODE with a step size of
% $latex \Delta t = .1$$. In other words,
% given $latex x_{k-1}$$, the value of the state at time $latex t_{k-1}$$,
% we model $latex x_k$$, the state value at time
% $latex t_k = t_{k-1} + \Delta t$$,
% by
% $latex \[
% \begin{array}{rcl}
%       x_{1,k} & = & x_{1,k-1} + x_{2,k-1} \Delta t + w_{1,k}
%       \\
%       x_{2,k} & = & x_{2,k-1} +
%       [ \mu ( 1 - x_{1,k-1}^2 ) x_{2,k-1} - x_{1,k-1} ] \Delta t
%       + w_{2,k}
% \end{array}
% \] $$
% where $latex w_k \sim \B{N}( 0 , Q_k )$$
% and the variance $latex Q_k = \R{diag} ( 0.01, 0.01 )$$.
% The routine $cref vanderpol_g.m$$ calculates the mean for $latex x_k$$
% given $latex x_{k-1}$$.
% Note that the initial state estimate is calculated by this routine as
% $latex \hat{x} = g_1 ( x_0 )$$ and has $latex \hat{x}_1 = \hat{x}_2 = 0$$
% and $latex Q_1 = \R{diag} ( 100. , 100. )$$.
%
% $head Measurements$$
% For this example, the measurements are noisy direct observations of position
% $latex \[
% \begin{array}{rcl}
%       z_{1,k} & = & x_{1,k} + v_{1,k}
% \end{array}
% \] $$
% where $latex v_k \sim \B{N}( 0 , R_k )$$ and
% the measurement variance $latex R_k = \R{diag} ( 0.01, 0.01 )$$.
% The routine $cref direct_h.m$$ calculates the mean for $latex z_k$$,
% given the value of $latex x_k$$; i.e. $latex x_{1,k}$$.
%
% $head Constraints$$
% There are no constraints for this example so we define
% $latex \[
%       f_k ( x_k ) = - 1
% \] $$
% The routine $cref no_f.m$$ is used for these calculations.
%
% $children%
%       example/nonlinear/vanderpol_sim.m%
%       example/nonlinear/vanderpol_g.m
% %$$
% $head Call Back Functions$$
% $table
% $rref vanderpol_sim$$
% $rref vanderpol_g.m$$
% $rref direct_h.m$$
% $rref no_f.m$$
% $tend
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = vanderpol_ok(draw_plot)
    if nargin < 1
        draw_plot = false;
    end
    % ---------------------------------------------------------------
    % Simulation problem parameters
    max_itr = 20;           % Maximum number of optimizer iterations
    epsilon = 1e-4;         % Convergence criteria
    N       = 41;           % Number of measurement points (N > 1)
    T       = 4.;           % Total time
    ell     = 1;            % Number of constraints (never active)
    n       = 2;            % Number of components in the state vector (n = 2)
    m       = 1;            % Number of measurements at each time point (m <= n)
    xi      = [ 0. , -5.]'; % True initial value for the state at time zero
    x_in    = zeros(n, N);  % Initial sequence where optimization begins
    initial = [ 0. , -5.]'; % Estimate of initial state (at time index one)
    mu      = 2.0;          % ODE nonlinearity parameter
    sigma_r = 1.;           % Standard deviation of measurement noise
    sigma_q = .1;           % Standard deviation of the transition noise
    sigma_e = 10;           % Standard deviation of initial state estimate
    randn('seed', 4321);    % Random number generator seed
    % -------------------------------------------------------------------
    % Level of tracing during optimization
    if draw_plot
        level = 1;
    else
        level = 0;
    end
    % spacing between time points
    dt = T / (N - 1);
    % -------------------------------------------------------------------
    %  information used by vanderpol_g
    params.vanderpol_g_initial  = initial;
    params.vanderpol_g_dt       = dt;
    params.vanderpol_g_mu       = mu;
    % -------------------------------------------------------------------
    % global information used by direct_h
    %global direct_h_info
    params.direct_h_index = 1;
    % ---------------------------------------------------------------
    % Rest of the information required by ckbs_nonlinear
    %
    % Step size for fourth order Runge-Kutta method in vanderpol_sim
    % is the same as time between measurement points
    step    = dt;
    % Simulate the true values for the state vector
    x_true  = vanderpol_sim(mu, xi, N, step);
    time    = linspace(0., T, N);
    % Simulate the measurement values
    z       = x_true(1:m,:) + sigma_r * randn(m, N);
    % Inverse covariance of the measurement and transition noise
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
    g_fun = @ (k,x) vanderpol_g(k,x,params);
    h_fun = @(k,x) direct_h(k,x,params);
    % ---------------------------------------------------------------
    % call the optimizer
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
        plot(time, x_out(1,:),  'g-');
        plot(time, z(1,:),      'ko');
        title('Position: blue=truth, green=estimate, o=data');
        hold off
        return
    end
end
% $$ $newlinech %$$
% $end
