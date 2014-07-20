% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin sine_wave_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       itr
%       dt
%       randn
%       cos
%       pos_vel
%       rinv
%       inv
%       qk
%       qinv
%       qinvk
%       complementarity
%       df
%       dg
%       dh
%       xk
%       uk
%       fk
%       gk
%       hk
%       sumsq
%       dk
%       clf
%       fid
%       fopen
%       dat
%       wt
%       fprintf
%       params
% $$
%
% $section ckbs_nonlinear: Example and Test of Tracking a Sine Wave$$
%
% $index ckbs_nonlinear, sine wave$$
% $index example, nonlinear$$
% $index test, nonlinear$$
%
%
% $head Syntax$$
% $syntax%[%ok%] = sine_wave_ok(%constraint%, %draw_plot%)%$$
%
% $head constraint$$
% is a character row vector with one of the following values:
% $table
% $bold value$$            $cnext $bold meaning$$                $rnext
% $code no_constraint$$    $cnext solve unconstrained problem    $rnext
% $code box_constraint$$   $cnext use box constraints            $rnext
% $code sine_constraint$$  $cnext use a sine wave constraint
% $tend
%
% $head draw_plot$$
% If $icode draw_plot$$ is false, no plot is drawn.
% If this argument is true, a plot is drawn showing the results.
%
% $head ok$$
% If the return value $italic ok$$ is true, the test passed,
% otherwise the test failed.
%
% $head Simulated State Vector$$
% $table
% $latex x1 (t) = 1$$               $cnext first component of velocity
% $rnext
% $latex x2 (t) = t$$               $cnext first component of position
% $rnext
% $latex x3 (t) = \cos(t)$$         $cnext second component of velocity
% $rnext
% $latex x4 (t) = \sin(t)$$         $cnext second component of position  
% $tend
%
% $head Transition Model$$
% We simulate imperfect knowledge of the dynamical system as follows
% $latex \[
% \begin{array}{rcl}
%       x_{1,k} & = & x_{1,k-1} + x_{2,k-1} \Delta t + w_{1,k}
%       \\
%       x_{2,k} & = & x_{2,k-1} + w_{2,k}
%       \\
%       x_{3,k} & = & x_{3,k-1} + x_{4,k-1} \Delta t + w_{3,k}
%       \\
%       x_{4,k} & = & x_{4,k-1} + w_{4,k}
% \end{array}
% \] $$
% where $latex \Delta t = 2 \pi / N$$ is the time between measurements,
% $latex N = 50$$ is the number of measurement times,
% $latex w_k \sim \B{N}( 0 , Q_k )$$,
% and the transition variance is given by
% $latex \[
% Q_k = \left( \begin{array}{cccc}
%       \Delta t       & \Delta t^2 / 2 & 0 & 0 \\
%       \Delta t^2 / 2 & \Delta t^3 / 3 & 0 & 0 \\
%       0              & 0              & \Delta t       & \Delta t^2 / 2 \\
%       0              & 0              & \Delta t^2 / 2 & \Delta t^3 / 3
% \end{array} \right)
% \] $$
% The routine $cref pos_vel_g.m$$ calculates the mean for $latex x_k$$
% given $latex x_{k-1}$$.
%
% $head Initial State Estimate$$
% The routine $cref pos_vel_g.m$$ also returns the initial state estimate 
% $latex \hat{x} = g_1 ( x_0 )$$ as the simulated state value at the 
% first measurement point; i.e.,
% $latex \[
%       \hat{x} = [ 1 , \Delta t , \cos ( \Delta t ) , \sin ( \Delta t ) ]^\R{T}
% \] $$
% The variance of  the initial state estimate is
% $latex Q_1 = 10^4 I_4 $$ where $latex I_4$$ is the four by four
% identity matrix.
%
% $head Measurement Model$$
% For this example, the measurements are noisy observations
% of the distance to two positions 
% $latex a = (0, -1.5)^\R{T}$$ and $latex b = ( 2 \pi , -1.5)$$; i.e.,
% $latex \[
% \begin{array}{rcl}
%       z_{1,k} & = & \sqrt{ ( x_{2,k} - a_1 )^2 + ( x_{4,k} - a_2 )^2 } + v_{1,k}
%       \\
%       z_{2,k} & = & \sqrt{ ( x_{2,k} - b_1 )^2 + ( x_{4,k} - b_2 )^2 } + v_{2,k}
% \end{array}
% \] $$
% where $latex v_k \sim \B{N} ( 0 , R_k )$$ and
% the measurement variance $latex R_k = \R{diag}( .25 , .25 )$$.
% The routine $cref distance_h.m$$ calculates the mean for $latex z_k$$,
% given the value of $latex x_k$$.
%
% $subhead Exception$$
% The simulated state sequence and the measurement model above are used
% to simulate the data with the following exception:
% if $latex \hat{k}$$ is the time index where $latex x_{4,k}$$ is maximal,
% $latex \[
% \begin{array}{rcl}
%       z_{1,k} & = & 
%       \sqrt{ ( x_{2,k} - a_1 )^2 + ( x_{4,k} + .5 - a_2 )^2 } 
%       \\
%       z_{2,k} & = & 
%       \sqrt{ ( x_{2,k} - b_1 )^2 + ( x_{4,k} + .5 - b_2 )^2 }
% \end{array}
% \] $$
% This increases the probability that the upper constraint will be
% active at the optimal solution.
%
% $head Constraint Model$$
% The are three possible constraints depending on the $icode constraint$$
% argument:
%
% $subhead no_constraint$$
% if $icode%constraint% == 'no_constraint'%$$, the function $cref no_f.m$$
% is used and the unconstrained problem is solved.
%
% $subhead box_constraint$$
% if $icode%constraint% == 'box_constraint'%$$, the function $cref box_f.m$$
% is used to represent the  following constraints:
% for $latex k = 1 , \ldots , N$$,
% $latex \[
%       -1 \leq x_{4,k} \leq +1
% \] $$
%
% $subhead sine_constraint$$
% if $icode%constraint% == 'sine_constraint'%$$, the function $cref sine_f.m$$
% is used to represent the  following constraints:
% for $latex k = 1 , \ldots , N$$,
% $latex \[
%       x_{4,k} \leq \sin( x_{2,k} ) + .1
% \] $$
%
%
% $head Call Back Functions$$
% $table
% $rref pos_vel_g.m$$
% $rref distance_h.m$$
% $rref no_f.m$$
% $rref box_f.m$$
% $rref sine_f.m$$
% $tend
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = sine_wave_ok(constraint, draw_plot)
    % --------------------------------------------------------
    % Simulation problem parameters
    max_itr    = 100;             % maximum number of iterations
    epsilon    = 1e-5;           % convergence criteria
    N          = 50;             % number of measurement time points
    dt         = 2 * pi / N;     % time between measurement points
    sigma_r    = .5;             % standard deviation of measurement noise
    sigma_q    = 1;              % square root of multiplier for transition variance
    sigma_e    = 100;            % standard deviation           
    delta_box  = 0.;             % distance from truth to box constraint
    delta_sine = .1;             % distance from truth to sine constraint
    h_min      = 0;              % minimum horizontal value in plots
    h_max      = 7;              % maximum horizontal value in plots
    v_min      = -1.5;           % minimum vertical value in plots
    v_max      = +1.5;           % maximum vertical value in plots 
    a          = [ 0    ; -1.5]; % position that range one measures from
    b          = [ 2*pi ; -1.5]; % position that range two measures to
    randn('seed', 1234);         % Random number generator seed
    % ---------------------------------------------------------
    % level of tracing during the optimization
    if draw_plot
            level   = 1;
    else
            level   = 0;
    end
    % ---------------------------------------------------------
    % simulate the true trajectory and measurement noise
    t        =  (1 : N) * dt;         % time
    x_true   = [ ones(1, N) ; t ; cos(t) ; sin(t) ];
    initial  = x_true(:,1);
    % ---------------------------------------------------------
    % information used by box_f
    params.box_f_lower = - 1 - delta_box;
    params.box_f_upper = + 1 + delta_box;
    params.box_f_index = 4;
    % ---------------------------------------------------------
    % information used by sine_f

    params.sine_f_index  = [ 2 ; 4 ];
    params.sine_f_offset = [ 0 ; delta_sine ];
    % ---------------------------------------------------------
    % information used by pos_vel_g
    params.pos_vel_g_dt       = dt;
    params.pos_vel_g_initial  = initial;
    g_fun     = @(k,x) pos_vel_g(k,x,params);

    % ---------------------------------------------------------
    %  information used by distance_h
   
    params.distance_h_index    = [ 2 ; 4 ];
    params.distance_h_position = [a , b ];
    h_fun = @(k,x) distance_h(k,x, params);
    % ---------------------------------------------------------
    % problem dimensions
    m     = 2; % number of measurements per time point
    n     = 4; % number of state vector components per time point
    ell   = 0; % number of constraint (reset to proper value below)
    % ---------------------------------------------------------
    % simulate the measurements
    [c,k_max] = max(x_true(4,:));    % index where second position component is max
    rinv      = zeros(m, m, N);      % inverse covariance of measurement noise
    z         = zeros(m, N);         % measurement values
    randn('state',sum(100*clock));   % random input
    for k = 1 : N
            x_k   = x_true(:, k);
            if k == k_max
                    % special case to increase chance that upper constraint is active 
                    x_k(2)   = x_k(2) + .5;
                    h_k      = h_fun(k, x_k);
                    z(:, k)  = h_k; 
            else
                    % simulation according to measurement model
                    h_k      = h_fun(k, x_k);
                    z(:, k)  = h_k + sigma_r * randn(2, 1); 
            end
            rinv(:,:, k) = eye(m) / (sigma_r * sigma_r);
    end

    %fid = fopen('z4.dat', 'wt');
    %fprintf(fid, '%6.6f %6.6f\n', z);
    z = load('z3.dat')';

    % ---------------------------------------------------------
    % inverse transition variance
    qk       = sigma_q * sigma_q * [ ...
            dt       , dt^2 / 2 , 0 , 0 
            dt^2 / 2 , dt^3 / 3 , 0 , 0 
            0        , 0        , dt       , dt^2 / 2
            0        , 0        , dt^2 / 2 , dt^3 / 3
    ];
    qinvk    = inv(qk);
    qinv     = zeros(n, n, N);
    for k = 2 : N
            qinv(:,:, k) = qinvk;
    end
    %
    % inverse initial estimate covariance
    qinv(:,:,1) = eye(n) * sigma_e * sigma_e; 
    % ----------------------------------------------------------
    % initialize x vector at position and velocity zero
    x_in = zeros(n, N);
    % ----------------------------------------------------------------------
   
    %
    is_no_f   = false;
    is_box_f  = false;
    is_sine_f = false;
    switch constraint
            case {'no_constraint'}
            ell       = 1;
            f_fun     = @ no_f;
            is_no_f   = true;

            case {'box_constraint'}
            ell       = 2;
            f_fun     = @(k,xk) box_f(k,xk,params);
            is_box_f  = true;

            case {'sine_constraint'}
            ell       = 1;
            f_fun     = @(k,xk) sine_f(k,xk,params);
            is_sine_f = true;

            otherwise
            constraint = constraint
            error('sine_wave_ok: invalid value for constraint')
    end
    % ----------------------------------------------------------------------
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
    ok         = true;
    ok         = ok & (size(info,1) <= max_itr);
    %
    % Set up affine problem (for change in x) corresponding to x_out.
    % Check constraint feasibility and complementarity.
    f_out      = zeros(ell, N);
    g_out      = zeros(n, N);
    h_out      = zeros(m, N);
    df_out     = zeros(ell, n, N);
    dg_out     = zeros(n, n, N);
    dh_out     = zeros(m, n, N);
    xk1        = zeros(n, 1);
    for k = 1 : N
            xk    = x_out(:, k);
            uk    = u_out(:, k);
            [fk, Fk]   = f_fun(k, xk);
            [gk, Gk]   = g_fun(k, xk1);
            [hk, Hk]   = h_fun(k, xk);
            %
            ok   = ok & all( fk <= epsilon );              % feasibility
            ok   = ok & all( abs(fk) .* uk <= epsilon );   % complementarity
            %
            f_out(:, k)    = fk;
            g_out(:, k)    = gk - xk;
            h_out(:, k)    = hk;
            df_out(:,:, k) = Fk;
            dg_out(:,:, k) = Gk;
            dh_out(:,:, k) = Hk;
            %
            xk1 = xk;
    end
    ok   = ok & all( all( u_out >= 0 ) );
    %
    % Compute gradient of unconstrained affine problem at zero
    % and check gradient of Lagrangian
    dx_out = zeros(n, N);
    d_out  = ckbs_sumsq_grad(dx_out, z, g_out, h_out, dg_out, dh_out, qinv, rinv);
    for k = 1 : N
            uk = u_out(:, k);
            Fk = df_out(:,:, k);
            dk = d_out(:, k);
            %
            ok = ok & (max ( abs( Fk' * uk + dk ) ) <= 1e2*epsilon);
    end
    %
    if draw_plot
            figure(1)
            clf
            hold on
            plot( x_true(2,:) , x_true(4,:) , 'b-' );
            plot( x_out(2,:)  , x_out(4,:)  , 'g-' );
            if is_no_f
                    title('blue=truth, green=estimate, no constraint')
            end
            if is_box_f
                    title('blue=truth, green=estimate, red=linear constraint')
                    plot( t , + ones(1, N) + delta_box , 'r-'); 
                    plot( t , - ones(1, N) - delta_box , 'r-'); 
            end
            if is_sine_f
                    title('blue=truth, green=estimate, red=nonlinear constraint')
                    plot( x_out(2,:) , sin( x_out(2,:) ) + delta_sine , 'r-'); 
            end
            axis( [ h_min, h_max, v_min, v_max] );
    return
end
% $$ $newlinech %$$
% $end
