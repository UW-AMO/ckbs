% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_nonlinear$$ $newlinech %$$
% $spell
%       itr
%       complementarity
%       obj
%       ckbs
%       qinv
%       rinv
%       feval
%       fk
%       gk
%       hk
%       xk
%       fu
%       aff
% $$
%
% $section The Nonlinear Constrained Kalman-Bucy Smoother$$
%
% $index ckbs_nonlinear$$
% $index constrained, Kalman-Bucy Smoother$$
% $index smoother, constrained Kalman-Bucy$$
% $index Kalman, constrained smoother$$
% $index Bucy, constrained smoother$$
%
% $head Syntax$$
% $codei/[/x_out/, /u_out/, /info/] = ckbs_nonlinear(/
%       f_fun/, /g_fun/, /h_fun/, ...
%       /max_itr/, /epsilon/, /x_in/, /z/, /qinv/, /rinv/, /level/)/$$
%
% $head Purpose$$
% This routine minimizes the
% Kalman-Bucy smoother residual sum of squares objective
% for a general nonlinear model functions
% (in the case where the model functions are affine, the routine
% $cref ckbs_affine$$ is more efficient).
%
% $head Notation$$
% The Kalman-Bucy smoother residual sum of squares is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = &
% \frac{1}{2}
% [ z_k - h_k ( x_k ) ]^\R{T} * R_k^{-1} * [ z_k - h_k ( x_k ) ]
% \\
% & + &
% \frac{1}{2}
% [ x_k - g_k ( x_{k-1} ) ]^\R{T} * Q_k^{-1} * [ x_k - g_k ( x_{k-1} ) ]
% \\
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
% Note that $latex g_1 ( x_0 )$$ is the initial state estimate
% and $latex Q_1$$ is the corresponding covariance.
%
% $head Problem$$
% Our constrained Kalman-Bucy smoother problem is
% $latex \[
% \begin{array}{rll}
% {\rm minimize}
%       & S( x_1 , \ldots , x_N )
%       & {\rm w.r.t.} \; x_1 \in \B{R}^n , \ldots , x_N \in \B{R}^n
% \\
% {\rm subject \; to}
%       & f_k(x_k) \leq 0
%       & {\rm for} \; k = 1 , \ldots , N
% \end{array}
% \] $$
%
% $head First Order Conditions$$
% A state sequence $latex ( x_1 , \ldots , x_N )$$ is considered a solution
% if there is a Lagrange multiplier sequence $latex ( u_1 , \ldots , u_N )$$
% such that the following conditions are satisfied for
% $latex i = 1 , \ldots , \ell$$ and $latex k = 1 , \ldots , N$$:
% $latex \[
% \begin{array}{rcl}
% f_k ( x_k )                                & \leq & \varepsilon  \\
% 0                                          & \leq & u_k          \\
% | u_k^\R{T} * \partial_k f_k ( x_k ) + \partial_k  S( x_1 , \ldots , x_N ) |
%                                            & \leq & \varepsilon \\
% | f_k ( x_k )_i | * ( u_k )_i              & \leq & \varepsilon
% \end{array}
% \] $$
% Here the notation
% $latex \partial_k$$ is used for the partial derivative of
% with respect to $latex x_k$$ and the notation
% $latex (u_k)_i$$ denotes the $th i$$ component of $latex u_k$$.
%
% $head f_fun$$
% The $code ckbs_nonlinear$$ argument $icode f_fun$$
% is a function that supports both of the
% following syntaxes
% $codei/
%       [/fk/] = feval(/f_fun/, /k/, /xk/)
%       [/fk/, /Fk/] = feval(/f_fun/, /k/, /xk/)
% /$$
%
% $subhead no_f$$
% In the case where there are no constraints,
% one can use $cref no_f.m$$ which is equivalent to
% $codep
%       [fk, Fk] = f_fun(k, xk)
%       n  = size(xk, 1);
%       fk = -1;
%       Fk =  zeros(1, n);
%       return
%       end
% $$
%
% $subhead k$$
% The $icode f_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$.
%
% $subhead xk$$
% The $icode f_fun$$ argument $icode xk$$ is an column vector with
% length $latex n$$. It specifies the state vector at time index $icode k$$
% $latex \[
%       xk = x_k
% \] $$.
%
% $subhead fk$$
% The $icode f_fun$$ result $icode fk$$ is an column vector with
% length $latex \ell$$ and
% $latex \[
%       fk = f_k ( x_k )
% \] $$
%
% $subhead Fk$$
% If the $icode f_fun$$ result $icode Fk$$ is present in the syntax,
% it is the $latex \ell \times n$$ matrix given by
% and
% $latex \[
%       Fk = \partial_k f_k ( x_k )
% \] $$
%
% $head g_fun$$
% The $code ckbs_nonlinear$$ argument $icode g_fun$$
% is a function that supports both of
% the following syntaxes
% $codei/
%       [/gk/] = feval(/g_fun/, /k/, /xk1/)
%       [/gk/, /Gk/] = feval(/g_fun/, /k/, /xk1/)
% /$$
%
% $subhead k$$
% The $icode g_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$.
% The case $latex k = 1$$ serves to specify the initial state estimate.
%
% $subhead xk1$$
% The $icode g_fun$$ argument $icode xk1$$ is an column vector with
% length $latex n$$. It specifies the state vector at time index $icode k$$
% $latex \[
%       xk1 = x_{k-1}
% \] $$.
% In the case $latex k = 1$$, the value of $icode xk1$$ does not matter.
%
% $subhead gk$$
% The $icode g_fun$$ result $icode gk$$ is an column vector with
% length $latex n$$ and
% $latex \[
%       gk = g_k ( x_{k-1} )
% \] $$
% In the case $latex k = 1$$, the value of $icode gk$$ is the initial state
% estimate at time index $icode k$$.
%
% $subhead Gk$$
% If the $icode g_fun$$ result $icode Gk$$ is present in the syntax,
% it is the $latex n \times n$$ matrix given by
% and
% $latex \[
%       Gk = \partial_{k-1} g_k ( x_{k-1} )
% \] $$
% In the case $latex k = 1$$, the value of $icode Gk$$ is the zero matrix;
% i.e., $latex g_k$$ does not depend on the value of $latex x_{k-1}$$.
%
% $head h_fun$$
% The $code ckbs_nonlinear$$ argument $icode h_fun$$
% is a function that supports both of the
% following syntaxes
% $codei/
%       [/hk/] = feval(/h_fun/, /k/, /xk/)
%       [/hk/, /Hk/] = feval(/h_fun/, /k/, /xk/)
% /$$
%
% $subhead k$$
% The $icode h_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$.
%
% $subhead xk$$
% The $icode h_fun$$ argument $icode xk$$ is an column vector with
% length $latex n$$. It specifies the state vector at time index $icode k$$
% $latex \[
%       xk = x_k
% \] $$.
%
% $subhead hk$$
% The $icode h_fun$$ result $icode hk$$ is an column vector with
% length $latex m$$ and
% $latex \[
%       hk = h_k ( x_k )
% \] $$
%
% $subhead Hk$$
% If the $icode h_fun$$ result $icode Hk$$ is present in the syntax,
% it is the $latex m \times n$$ matrix given by
% and
% $latex \[
%       Hk = \partial_k h_k ( x_k )
% \] $$
%
% $head max_itr$$
% The integer scalar $icode max_itr$$ specifies the maximum number of
% iterations of the algorithm to execute. It must be greater than or
% equal to zero. Note that if it is zero, the first row of the
% $cref/info/ckbs_nonlinear/info/$$ return value will still be computed.
% This can be useful for deciding what is a good value for the argument
% $cref/epsilon/ckbs_nonlinear/epsilon/$$.
%
% $head epsilon$$
% The $code ckbs_nonlinear$$ argument $icode epsilon$$ is a positive scalar.
% It specifies the convergence
% criteria value; i.e.,
% $latex \[
%       \varepsilon = epsilon
% \] $$
% In some instances, $code ckbs_nonlinear$$ will return after printing
% a warning without convergence; see $cref/info/ckbs_nonlinear/info/$$.
%
% $head x_in$$
% The $code ckbs_nonlinear$$ argument $icode x_in$$
% is a two dimensional array with size $latex n \times N$$.
% It specifies a sequence of state values; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x\_in (:, k) = x_k
% \] $$
% The closer the initial state sequence is to the solution
% the faster, and more likely, the $code ckbs_nonlinear$$ will converge.
% The initial state sequence need not be feasible; i.e.
% it is not necessary that
% $latex \[
%       f_k ( x_k ) \leq 0
% \] $$
% for all $latex k$$.
%
% $head z$$
% The $code ckbs_nonlinear$$ argument $icode z$$ is a two dimensional array
% with size $latex m \times N$$.
% It specifies the sequence of measurement vectors; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       z(:, k) = z_k
% \]$$
%
% $head qinv$$
% The $code ckbs_nonlinear$$ argument $icode qinv$$
% is a three dimensional array
% with size $latex n \times n \times N$$.
% It specifies the inverse of the variance of the measurement noise; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       qinv(:,:, k) = Q_k^{-1}
% \]$$
% In the case $latex k = 1$$, the value of $latex Q_k$$ is the variance
% of the initial state estimate (see $cref/g_fun/ckbs_nonlinear/g_fun/$$.
%
% $head rinv$$
% The $code ckbs_nonlinear$$ argument $icode rinv$$
% is a three dimensional array,
% with size $latex m \times m \times N$$.
% it specifies the inverse of the variance of the transition noise; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       rinv(:,:, k) = R_k^{-1}
% \]$$
% It is ok to signify a missing data value by setting the corresponding
% row and column of $icode rinv$$ to zero. This also enables use
% to interpolate the state vector $latex x_k$$ to points where
% there are no measurements.
%
% $head level$$
% The integer scalar $icode level$$ specifies the level of tracing to do.
% If $icode%level% == 0%$$, no tracing is done.
% If $icode%level% >= 1%$$, the row index $icode q$$
% and the correspond row of $icode info$$
% are printed as soon as soon as they are computed.
% If $icode%level% >= 2%$$, a check of the derivative calculations
% is printed before the calculations start. In this case, control will
% return to the keyboard so that you can print values, continue,
% or abort the calculation.
%
% $head x_out$$
% The result $icode x_out$$
% is a two dimensional array with size $latex n \times N$$.
% $latex ( x_1 , \ldots , x_N )$$.
% It contains an approximation for the optimal sequence; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x\_out(:, k) = x_k
% \]$$
%
% $head u_out$$
% The result $icode u_out$$
% is a two dimensional array with size % $latex \ell \times N$$.
% It contains the Lagrange multiplier sequence
% corresponding to $icode x_out$$; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       u\_out(:, k) = u_k
% \]$$
% The pair $icode x_out$$, $icode u_out$$ satisfy the
% $cref/first order conditions/ckbs_nonlinear/First Order Conditions/$$.
%
% $head info$$
% The result $icode info$$ is a matrix with each row corresponding
% to an iteration of the algorithm.
% Note that $code ckbs_nonlinear$$ has satisfied the convergence condition if
% and only if
% $codei%
%       all( %info%(end, 1:3) <= %epsilon% )
% %$$
% $pre
%
% $$
% We use $latex \{ x_k^q \}$$ to denote the state vector sequence at the
% at the end of iteration $latex q-1$$ and
% $latex \{ u_k^q \}$$ for the dual vector sequence ($latex u_k^q \geq 0$$).
%
% $subhead max_f$$
% The value
% $codei%
%       %max_f%(%q%) = %info%(%q%, 1)
% %$$
% is a bound on the constraint functions
% at the end of iteration $icode q-1$$. To be specific
% $latex \[
%        max\_f(q) = \max_{i, k} [ f_k ( x_k^q )_i ]
% \] $$
% where the maximum is with respect to
% $latex i = 1 , \ldots , \ell$$ and
% $latex k = 1 , \ldots , N$$.
% Convergence requires this value to be less than or equal $icode epsilon$$.
%
% $subhead max_grad$$
% The value
% $codei%
%       %max_grad%(%q)% = %info%(%q%, 2)
% %$$
% is a bound on the partial of the
% Lagrangian with respect to $latex x$$
% at the end of iteration $icode q-1$$. To be specific
% $latex \[
% max\_grad(q) = \max_{j, k} \left[ \left|
% (u_k^q)^\R{T} * \partial_k f_k(x_k^q) + \partial_k  S(x_1^q, \ldots ,x_N^q)
% \right| \right]_j
% \] $$
% where the maximum is with respect to
% $latex j = 1 , \ldots , n$$ and
% $latex k = 1 , \ldots , N$$.
% This value should converge to zero.
%
% $subhead max_fu$$
% The value
% $codei%
%       %max_fu% = %info%(%q%, 3)
% %$$
% is a bound on
% the complementarity conditions
% at the end of iteration $icode q-1$$. To be specific
% $latex \[
% max\_fu(q) = \max_{i, k} [ | f_k ( x_k )_i | * ( u_k )_i ]
% \] $$
% where the maximum is with respect to
% $latex i = 1 , \ldots , \ell$$ and
% $latex k = 1 , \ldots , N$$.
% This value should converge to be less than or equal zero.
%
% $subhead obj_cur$$
% The value
% $codei%
%       %obj_cur%(%q%) = %info%(%q%, 4)
% %$$
% is the value of the objective function
% at the end of $latex q-1$$; i.e.,
% $latex \[
%       obj\_cur(q) = S( x_1^q , \ldots , x_k^q )
% \] $$
% (Note that $latex info(1, 4)%$$ is the value of the objective
% corresponding to $icode x_in$$.
%
% $subhead lambda$$
% The value
% $codei%
%       %lambda%(%q%) = %info%(%q%, 5)
% %$$
% is the line search step size used during
% iteration $latex q-1$$ of $code ckbs_nonlinear$$.
% If the problem is nearly affine (if the affine approximate is accurate)
% this will be one.
% Much small step sizes indicate highly non-affine problems
% (with the exception that $icode%info%(1, 5)%$$ is always zero).
%
% $subhead lam_aff$$
% The value
% $codei%
%       %lam_aff%(%q%) = %info%(%q%, 6)
% %$$
% is the minimum line search step size
% while solving the affine sub-problem that corresponds to
% iteration $latex q-1$$ of $code ckbs_nonlinear$$.
% If the sub-problem solution is working well,
% this will be one.
% Much small step sizes indicate the sub-problem solution is not working well.
%
% $subhead Penalty Parameter$$
% The exact penalty function
% $latex \[
%       S( x_1 , \ldots , x_k )
%       + \alpha \sum_{k=1}^N \sum_{i=1}^\ell \max [ f_k ( x_k )_i , 0 ]
% \] $$
% is used to determine the line search step size.
% The value $icode%info%(%q%, 7)%$$ is the penalty parameter
% $latex \alpha$$ during iteration $latex q-1$$ of $code ckbs_nonlinear$$.
%
% $children%
%       example/nonlinear/get_started_ok.m%
%       example/nonlinear/sine_wave_ok.m%
%       example/nonlinear/vanderpol_ok.m%
%       example/nonlinear/vanderpol_experiment_simple.m%
%       example/nonlinear/utility.omh
% %$$
%
% $head Example$$
%
% $subhead Simple$$
% The file $cref get_started_ok.m$$ contains a simple example
% and test of $code ckbs_nonlinear$$.
% It returns true if $code ckbs_nonlinear$$ passes the tests
% and false otherwise.
%
% $subhead Unconstrained$$
% The option $icode%constraint% = 'no_constraint'%$$ to the routine
% $cref sine_wave_ok.m$$ runs an unconstrained example / test.
%
% $subhead Box Constraints$$
% The option $icode%constraint% = 'box_constraint'%$$ to the routine
% $cref sine_wave_ok.m$$ runs a example / test
% with upper and lower constraints.
%
% $subhead Nonlinear Constraints$$
% The option $icode%constraint% = 'sine_constraint'%$$ to the routine
% $cref sine_wave_ok.m$$ runs a example / test with a nonlinear constraint.
%
% $subhead Van Der Pol Oscillator$$
% The file $cref vanderpol_ok.m$$ contains an example
% and test of $code ckbs_nonlinear$$ estimating the position of a
% Van Der Pol oscillator (an oscillator with non-linear dynamics).
% It returns true if $code ckbs_nonlinear$$ passes the tests
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------
function [x_out, u_out, info] = ckbs_nonlinear(f_fun, g_fun, h_fun, ...
    max_itr, epsilon, x_in, z, qinv, rinv, level)
%

    %epsilon = epsilon/20;

    if nargin ~= 10
        error('ckbs_nonlinear: improper number of input arguments');
    end
    if nargout ~= 3
        error('ckbs_nonlinear: improper number of return values');
    end
    %
    % need n to evaluate gk
    n         = size(x_in, 1);
    %
    % a positve value smaller than epsilon
    delta = 1e-1 * epsilon;
    %
    [fk, Fk] = feval(f_fun, 1, x_in(:,1));
    [gk, Gk] = feval(g_fun, 1, zeros(n, 1));
    [hk, Hk] = feval(h_fun, 1, x_in(:,1));
    %
    if level == 2
        % initial state esitmate is used for checking derivatives
        x_initial_estimate = gk;
    end
    %
    % get other sizes of problem
    N     = size(x_in, 2);
    m     = size(z, 1);
    ell   = size(fk, 1);
    %
    % check for nan or infinity
    if ~ all( isfinite( x_in(:) ) ) | ...
            ~ all( isfinite( qinv(:) ) ) | ...
            ~ all( isfinite( rinv(:) ) )
        error('ckbs_nonlinear: argument is not finite valued');
    end
    if ~ all( isfinite( fk(:) ) ) | ~ all( isfinite( Fk(:) ) ) | ...
            ~ all( isfinite( gk(:) ) ) | ~ all( isfinite( Gk(:) ) ) | ...
            ~ all( isfinite( hk(:) ) ) | ~ all( isfinite( Hk(:) ) )
        error('ckbs_nonlinear: a callback function retrun is not finite valued');
    end
    %
    % check n in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if n ~= size(x_in,1) | n ~= size(qinv,1) | n ~= size(qinv,2) | ...
            n ~= size(Fk,2) | n ~= size(gk,1) | n ~= size(Gk,1) | n ~= size(Gk,2) | ...
            n ~= size(Hk,2)
        size(x_in,1)
        size(qinv,1)
        size(qinv,2)
        size(Fk,2)
        size(gk,1)
        size(Gk,1)
        size(Gk,2)
        size(Hk,2)
        error('ckbs_nonlinear: argument sizes with value n do not agree');
    end
    % check m in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if m ~= size(z,1) | m ~= size(rinv,1) | m ~= size(rinv,2) | ...
            m ~= size(hk,1) | m ~= size(Hk,1)
        size(z,1)
        size(rinv,1)
        size(rinv,2)
        size(hk,1)
        size(Hk,1)
        error('ckbs_nonlinear:  argument sizes with value m do not agree');
    end
    % check ell in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if ell ~= size(fk,1) | ell ~= size(Fk,1)
        size(fk,1)
        size(Fk,1)
        error('ckbs_nonlinear:  argument sizes with value ell do not agree');
    end
    % check N in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if N ~= size(x_in,2) | N ~= size(z,2) | N ~= size(qinv,3) | N ~= size(rinv,3)
        size(x_in,2)
        size(z,2)
        size(qinv,3)
        size(rinv,3)
        error('ckbs_nonlinear:  argument sizes with value N do not agree');
    end
    %
    % compare derivatives computed by f_fun h_fun and g_fun
    % with central difference approximations
    if level >= 2
        step = 1e-5;
        [f1, F1] = feval(f_fun, 1, x_initial_estimate);
        [g1, G1] = feval(g_fun, 2, x_initial_estimate);
        [h1, H1] = feval(h_fun, 1, x_initial_estimate);
        str      = sprintf('%12s%12s%12s%12s', ...
            'name', 'analytic', 'numerical', 'difference' ...
            );
        disp(str);
        for j = 1 : n
            x1     = x_initial_estimate(:,1);
            x1(j)  = x_initial_estimate(j,1) + step;
            fp     = feval(f_fun, 1, x1);
            gp     = feval(g_fun, 2, x1);
            hp     = feval(h_fun, 1, x1);
            %
            x1(j)  = x_initial_estimate(j,1) - step;
            fm     = feval(f_fun, 1, x1);
            hm     = feval(h_fun, 1, x1);
            gm     = feval(g_fun, 2, x1);
            %
            for i = 1 : ell
                name       = sprintf('df_%d/dx_%d', i, j);
                analytic   = F1(i, j);
                numerical  = (fp(i) - fm(i)) / (2*step);
                difference = analytic - numerical;
                str        = sprintf('%12s%12.3g%12.3g%12.3g', ...
                    name, analytic, numerical, difference ...
                    );
                disp(str);
            end
            for i = 1 : n
                name       = sprintf('dg_%d/dx_%d', i, j);
                analytic   = G1(i, j);
                numerical  = (gp(i) - gm(i)) / (2*step);
                difference = analytic - numerical;
                str        = sprintf('%12s%12.3g%12.3g%12.3g', ...
                    name, analytic, numerical, difference ...
                    );
                disp(str);
            end
            for i = 1 : m
                name       = sprintf('dh_%d/dx_%d', i, j);
                analytic   = H1(i, j);
                numerical  = (hp(i) - hm(i)) / (2*step);
                difference = analytic - numerical;
                str        = sprintf('%12s%12.3g%12.3g%12.3g', ...
                    name, analytic, numerical, difference ...
                    );
                disp(str);
            end
        end
        disp('ckbs_nonlinear: end derivative check')
        keyboard
    end
    %
    % a useful constant
    zero_n     = zeros(n, 1);
    zero_n_N   = zeros(n, N);
    %
    % dimension current arrays
    f_cur  = zeros(ell, N);
    g_cur  = zeros(n, N);
    h_cur  = zeros(m,  N);
    %
    df_cur = zeros(ell, n, N);
    dg_cur = zeros(n, n, N);
    dh_cur = zeros(m, n, N);
    %
    %  initialize some values
    x_cur      = x_in;        % best state sequence so far
    info.iters       = zeros(0, 7); % return convergence tracking information
    itr        = 0;           % iteration counter
    lambda     = 0;           % linear search step size
    max_affine = 30;          % maximum number of iterations in affine sub-problem
    alpha      = 0;           % initial exact penalty function parameter
    %
    % initialize current value
    xk1      = zero_n;
    dist_cur = 0;
    for k  = 1 : N
        xk                           = x_cur(:, k);
        [fk, Fk]                     = feval(f_fun, k, xk);
        [gk, Gk]                     = feval(g_fun, k, xk1);
        [hk, Hk]                     = feval(h_fun, k, xk);
        %
        dist_cur                     = dist_cur + sum( max(fk, 0 ) );
        %
        f_cur(:, k)                  = fk;
        g_cur(:, k)                  = gk - xk;
        h_cur(:, k)                  = hk;
        %
        df_cur(:,:, k)               = Fk;
        dg_cur(:,:, k)               = Gk;
        dh_cur(:,:, k)               = Hk;
        %
        xk1                          = xk;
    end
    %
    % algorithm iteration loop
    for itr = 0 : max_itr
        %
        % objective fucntions corresponding to x_cur
        obj_cur = ckbs_sumsq_obj(zero_n_N, z, ...
            g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
        %
        % gradient of objective corresponding to x_cur
        grad_cur = ckbs_sumsq_grad(zero_n_N, z, ...
            g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
        %
        % affine approximate subproblem
        [dx_cur, u_new, affine_info] = ckbs_affine(max_affine, delta, z, ...
            f_cur, g_cur, h_cur, df_cur, dg_cur, dh_cur, qinv, rinv);
        if( any( any( u_new < 0 ) ) ) ...
                error('ckbs_nonlinear:  u_new has a negative element');
        end
        if( any( affine_info.iters(end, 1:3) > delta ) )
            'warning: ckbs_nonlinear: affine sub-problem did not converge';
        end
        %
        % affine approximation to objective corresponding to
        % x_new = x_cur + dx_cur
        obj_aff = ckbs_sumsq_obj(dx_cur, z, ...
            g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
        %
        % compute the convergence information for the current iterate
        dist_aff = 0;
        for k = 1 : N
            uk          = u_new(:, k);
            xk          = x_cur(:, k);
            Fk          = df_cur(:,:, k);
            fk          = f_cur(:, k);
            uk_fk       = uk .* fk;
            %
            dk          = grad_cur(:, k);
            Fk_uk_dk    = Fk' * uk + dk;
            if( k == 1 )
                max_f    = max(fk);
                max_grad = max(abs(Fk_uk_dk));
                max_fu   = max(abs(uk_fk));
            end
            max_f       = max(max_f,     max(fk)             );
            max_grad    = max(max_grad,  max(abs(Fk_uk_dk))  );
            max_fu      = max(max_fu,    max(abs(uk_fk))     );
            %
            fk_new      = f_cur(:, k) + Fk * dx_cur(:, k);
            uk_fk_new   = uk .* fk_new;
            %
            % update distance of affine solution from constraint set
            dist_aff    = dist_aff + sum( max( fk_new , 0) );
        end
        % minimum step size used by affine sub-problem
        if( size(affine_info.iters, 1) > 1 )
            lam_aff  = min( affine_info.iters(2:end,4) );
        else
            lam_aff  = 1;
        end
        %
        % add informaiton for this iteration
        info_itr = [max_f, max_grad, max_fu, obj_cur, lambda, lam_aff, alpha];
        if level >= 1
            msg = [ ...
                'q = ', num2str(itr), ...
                ', max_f, max_grad, max_fu', ...
                ',  obj_cur, lambda, lam_aff, alpha =' ...
                ];
            disp(msg);
            disp( num2str(info_itr, '%10.3g') );
        end
        info.iters     = [ info.iters ; info_itr ];
        %
        % check for convergence
        converge = ( (max_f <= epsilon) ...
            & (max_grad <= epsilon) ...
            & (max_fu <= epsilon) );
        %
        % check for done
        if (itr == max_itr) | converge
            x_out   = x_cur;
            u_out   = u_new;
            info.D = affine_info.D;
            info.A = affine_info.A;
            return;
        end
        % ---------------------------------------------------------------
        % Line Search
        % ---------------------------------------------------------------
        if( dist_aff > epsilon )
            warning('ckbs_nonlinear: affine constraint not satisfied');
            dist_aff = dist_aff
            x_out   = x_cur;
            u_out   = u_new;
            return;
        end
        %
        % update alpha
        if alpha == 0
            alpha = max( u_new(:) ) / 10;
        end
        obj_dir_der   = grad_cur(:)' * dx_cur(:);
        quadratic_aff = obj_aff - obj_cur - obj_dir_der;
        gamma         = 2 * quadratic_aff + obj_dir_der;
        if (dist_cur > 0) & (gamma > alpha * dist_cur)
            ratio      = gamma / dist_cur;
            alpha      = max( ratio , 2 * alpha);
        end
        rate        = obj_dir_der - alpha * max(dist_cur - dist_aff, 0);
        total_cur   = obj_cur + alpha * dist_cur;
        %
        % line search
        %
        ok        = 0;
        kount     = 0;
        max_kount = 30;
        lambda    = 2.;
        while (~ ok ) & (kount < max_kount)
            kount  = kount + 1;
            lambda = lambda / 2;
            %
            % argument corresponding to this step
            x_lam = x_cur + lambda .* dx_cur;
            %
            % calculate  f, g, h, df, dg, and dh corresponding to x_lam
            % also use f_cur, g_cur, h_cur, df_cur, dg_cur, dh_cur with
            % to store f_lam, g_lam, h_lam, df_lam, dg_lam, dh_lam
            xk1      = zero_n;
            dist_lam = 0;
            for k      = 1 : N
                xk                 = x_lam(:,k);
                %
                [fk, Fk]           = feval(f_fun, k, xk);
                [gk, Gk]           = feval(g_fun, k, xk1);
                [hk, Hk]           = feval(h_fun, k, xk);
                %
                dist_lam           = dist_lam + sum( max(fk, 0) );
                %
                f_cur(:, k)        = fk;
                g_cur(:, k)        = gk - xk;
                h_cur(:,k)         = hk;
                %
                df_cur(:,:, k)     = Fk;
                dg_cur(:,:, k)     = Gk;
                dh_cur(:,:, k)     = Hk;
                %
                xk1                = xk;
            end
            %
            % objective corresponding to x_lam
            obj_lam = ckbs_sumsq_obj(zero_n_N, z, ...
                g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
            %
            % check if we need to increase penalty parameter
            total_lam  = obj_lam + alpha * dist_lam;
            %
            diff_lam = total_lam - total_cur;
            ok       = (diff_lam <= (lambda / 10) * rate);
        end
        if ~ok
            rate    = rate
            warning('ckbs_nonlinear: line search failed');
            x_out   = x_cur;
            u_out   = u_new;
            return;
        end
        %
        % advance for next iteration
        x_cur    = x_lam;
        dist_cur = dist_lam;
    end
    error('ckbs_nonlinear: program error');
    return
end
