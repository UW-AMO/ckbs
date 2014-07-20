% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_L1_nonlinear$$ $newlinech %$$
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
%       linearized
%       subproblem
%       iteratively
%       relinearizing
%       linearization
%       optimality
%       pairwise
%       kuhn
% $$
%
% $section The Nonlinear Constrained Kalman-Bucy Smoother$$
%
% $index ckbs_L1_nonlinear$$
% $index robust L1, Kalman Smoother$$
% $index smoother, robust L1 Kalman$$
% $index Kalman, robust L1 smoother$$
%
% $head Syntax$$
% $codei/[/xOut/, /rOut/, /sOut/, /pPlusOut/, /pMinusOut/, /info/] = ckbs_L1_nonlinear(/
%       f_fun/, /g_fun/, /h_fun/, ...
%       /max_itr/, /epsilon/, /x_in/, /z/, /qinv/, /rinv/, /level/)/$$
%
% $head Purpose$$
% This routine minimizes the robust L1
% Kalman smoother objective
% for general nonlinear process and measurement models
% (in the case where the model functions are affine, the routine
% $cref ckbs_L1_affine$$ is more efficient).
%
% $head Notation$$
% The state sequence $latex x $$ is defined by $latex x  = x_1, \ldots, x_N
% $$,
% with each $latex x_i \in \B{R}^n $$.
% The robust L1 Kalman-Bucy smoother objective is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x ) & = & S^Q ( x ) + S^R ( x )
% \\
%  S^Q ( x ) & = & \frac{1}{2}  \sum_{k=1}^N
% ( x_k - g_k ( x_{k-1} ) )^\R{T} * Q_k^{-1} * ( x_k - g_k ( x_{k-1} ) )
% \\
% S^R( x ) & = &  \sqrt{2} \sum_{k=1}^N  \|R_k^{-1/2}( z_k - h_k ( x_k  ) )\|_1
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
% Note that $latex g_1$$ is the initial state estimate
% and $latex Q_1$$ is the corresponding covariance.
%
% $head Problem$$
% The robust L1 Kalman smoother problem is
% $latex \[
% {\rm minimize}\;
%        S( x ) \;
%        {\rm w.r.t.} \; x
% \] $$
%
% $head Optimality Conditions $$
% Define $latex d = d_1, \ldots, d_N $$ with each $latex d_i \in \B{R}^n$$.
% Define
% $latex \[
% \begin{array}{rcl}
% \tilde{S} ( x ; d) & = & \tilde{S}^Q ( x ; d ) + \tilde{S}^R ( x
% ; d )
% \\
%  \tilde{S}^Q ( x ; d ) & = & \frac{1}{2}  \sum_{k=1}^N
% ( x_k - g_k ( x_{k-1} ) - \partial_k g_k ( x_{k-1} ) d_k ) ^\R{T} * Q_k^{-1} * ( x_k - g_k (
% x_{k-1} ) - \partial_k g_k ( x_{k-1} ) d_k )
% \\
% \tilde{S}^R( x ; d ) & = &  \sqrt{2}   \sum_{k=1}^N
% \|R_k^{-1/2}( z_k - h_k ( x_k  ) - \partial_k h_k ( x_k) d_k )\|_1
% \end{array}
% \] $$
% A state sequence $latex x  $$ is considered a
% solution if
% $latex \[
% S( x ) - \min_d\; \tilde{S} ( x ; d  )  <  \varepsilon \;.
% \]$$
% This follows directly from the convex composite structure of
% $latex S $$.
% $head f_fun$$
% Constraints are current not supported in the robust L1 smoother.
% The syntax was left unchanged to match the constrained smoother,
% and because in future versions constraints may be supported.
%
% $head g_fun$$
% The $code ckbs_L1_nonlinear$$ argument $icode g_fun$$
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
% The $code ckbs_L1_nonlinear$$ argument $icode h_fun$$
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
% $cref/info/ckbs_L1_nonlinear/info/$$ return value will still be computed.
% This can be useful for deciding what is a good value for the argument
% $cref/epsilon/ckbs_L1_nonlinear/epsilon/$$.
%
% $head epsilon$$
% The $code ckbs_L1_nonlinear$$ argument $icode epsilon$$ is a positive scalar.
% It specifies the convergence
% criteria value; i.e.,
% $latex \[
%       \varepsilon = epsilon
% \] $$
%
% $head x_in$$
% The $code ckbs_L1_nonlinear$$ argument $icode x_in$$
% is a two dimensional array with size $latex n \times N$$.
% It specifies a sequence of state values; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x\_in (:, k) = x_k
% \] $$
% The closer the initial state sequence is to the solution
% the faster, and more likely, the $code ckbs_L1_nonlinear$$ will converge.
%
% $head z$$
% The $code ckbs_L1_nonlinear$$ argument $icode z$$ is a two dimensional array
% with size $latex m \times N$$.
% It specifies the sequence of measurement vectors; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       z(:, k) = z_k
% \]$$
%
% $head qinv$$
% The $code ckbs_L1_nonlinear$$ argument $icode qinv$$
% is a three dimensional array
% with size $latex n \times n \times N$$.
% It specifies the inverse of the variance of the measurement noise; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       qinv(:,:, k) = Q_k^{-1}
% \]$$
% In the case $latex k = 1$$, the value of $latex Q_k$$ is the variance
% of the initial state estimate (see $cref/g_fun/ckbs_L1_nonlinear/g_fun/$$.
%
% $head rinv$$
% The $code ckbs_L1_nonlinear$$ argument $icode rinv$$
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
% $head rOut$$
% The result $icode rOut$$ contains the structure $icode rOut$$
% reported by $cref ckbs_L1_affine$$ for the last affine
% subproblem.
%
% $head sOut$$
% The result $icode sOut$$ contains the structure $icode sOut$$
% reported by $cref ckbs_L1_affine$$ for the last affine
% subproblem.
%
% $head pPlusOut$$
% The result $icode pPlusOut$$ contains the structure $icode pPlusOut$$
% reported by $cref ckbs_L1_affine$$ for the last affine
% subproblem.
%
% $head pMinusOut$$
% The result $icode pMinusOut$$ contains the structure $icode pMinusOut$$
% reported by $cref ckbs_L1_affine$$ for the last affine
% subproblem.
%
% $head info$$
% The result $icode info$$ is a matrix with each row corresponding
% to a solution of $cref ckbs_L1_affine$$ on an affine subproblem.
% There are five numbers output per iteration:
% the infinity norm of the Kuhn-Tucker conditions,
% the one-norm of the Kuhn-Tucker conditions, the affine improvement,
% $latex \mu$$, and the number of line search iterations.
% See $cref ckbs_L1_affine$$ for the definitions of these quantities.
% Note that $code ckbs_L1_affine$$ has satisfied the convergence
% condition if the affine improvement is close to zero, that is,
% $codei%
%        %info%(end, 3) <= %epsilon%
% %$$
% $pre
%
% $$
%
%
% $children%
%       example/nonlinear/L1_nonlinear_ok.m
% %$$
%
% $head Example$$
%
% The file $cref L1_nonlinear_ok.m$$ contains an example
% and test of $code ckbs_L1_nonlinear$$ estimating the position of a
% Van Der Pol oscillator (an oscillator with non-linear dynamics).
% The code runs both $cref ckbs_nonlinear$$ and $cref
% ckbs_L1_nonlinear$$ on the example, in cases where there are
% outliers in the data, and can be used to make a plot.
% It returns true if $code ckbs_L1_nonlinear$$ converged on the problem,
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------
function [x_out, r_out, s_out, p_plus_out, p_minus_out, info] = ckbs_L1_nonlinear(f_fun, g_fun, h_fun, ...
    max_itr, epsilon, x_in, z, qinv, rinv, level)
%
    if nargin ~= 10
        error('ckbs_L1_nonlinear: improper number of input arguments');
    end
    if nargout ~= 6
        error('ckbs_L1_nonlinear: improper number of return values');
    end

    % need n to evaluate gk
    n         = size(x_in, 1);
    %
    % a positve value smaller than epsilon
    %delta = 1e-2 * epsilon;
    delta = inf;  % only need to reduce duality gap.


    %
    [gk, Gk] = feval(g_fun, 1, zeros(n, 1));
    [hk, Hk] = feval(h_fun, 1, x_in(:,1));
    %

    if level == 2
        % initial state esitmate is used for checking derivatives
        x_initial_estimate = gk;
    end


    % get other sizes of problem
    N     = size(x_in, 2);
    m     = size(z, 1);
    %

    % check for nan or infinity
    if ~ all( isfinite( x_in(:) ) ) | ...
            ~ all( isfinite( qinv(:) ) ) | ...
            ~ all( isfinite( rinv(:) ) )
        error('ckbs_L1_nonlinear: argument is not finite valued');
    end
    if ~ all( isfinite( gk(:) ) ) | ~ all( isfinite( Gk(:) ) ) | ...
            ~ all( isfinite( hk(:) ) ) | ~ all( isfinite( Hk(:) ) )
        error('ckbs_L1_nonlinear: a callback function retrun is not finite valued');
    end


    % check n in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if n ~= size(x_in,1) || n ~= size(qinv,1) || n ~= size(qinv,2) || ...
            n ~= size(gk,1) || n ~= size(Gk,1) || n ~= size(Gk,2) || ...
            n ~= size(Hk,2)
        size(x_in,1)
        size(qinv,1)
        size(qinv,2)
        size(gk,1)
        size(Gk,1)
        size(Gk,2)
        size(Hk,2)
        error('ckbs_L1_nonlinear: argument sizes with value n do not agree');
    end
    % check m in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if m ~= size(z,1) || m ~= size(rinv,1) || m ~= size(rinv,2) || ...
            m ~= size(hk,1) || m ~= size(Hk,1)
        size(z,1)
        size(rinv,1)
        size(rinv,2)
        size(hk,1)
        size(Hk,1)
        error('ckbs_L1_nonlinear:  argument sizes with value m do not agree');
    end
    % check N in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
    if N ~= size(x_in,2) || N ~= size(z,2) || N ~= size(qinv,3) || N ~= size(rinv,3)
        size(x_in,2)
        size(z,2)
        size(qinv,3)
        size(rinv,3)
        error('ckbs_L1_nonlinear:  argument sizes with value N do not agree');
    end
    %
    % compare derivatives computed by f_fun h_fun and g_fun
    % with central difference approximations
    if level >= 2
        step = 1e-5;
        [g1, G1] = feval(g_fun, 2, x_in(:,1));
        [h1, H1] = feval(h_fun, 1, x_in(:,1));
        disp('analytic , numerical, difference');
        for j = 1 : n
            x1     = x_in(:,1);
            x1(j)  = x_in(j,1) + step;
            gp     = feval(g_fun, 2, x1);
            hp     = feval(h_fun, 1, x1);
            %
            x1(j)  = x_in(j,1) - step;
            hm     = feval(h_fun, 1, x1);
            gm     = feval(g_fun, 2, x1);
            %
            analytic   = [ G1(:, j) ; H1(:, j) ];
            numerical  = [ (gp-gm) ; (hp-hm) ] / (2*step);
            difference = analytic - numerical;
            disp( [analytic , numerical , difference ] );
        end
        keyboard('ckbs_L1_nonlinear: end derivative check')
    end
    %
    % a useful constant
    zero_n     = zeros(n, 1);
    zero_n_N   = zeros(n, N);
    %
    % dimension of array
    %
    g_cur  = zeros(n, N);
    h_cur  = zeros(m,  N);
    dg_cur = zeros(n, n, N);
    dh_cur = zeros(m, n, N);
    %
    %  initialize some values
    x_cur      = x_in;        % current state estimate
    info       = zeros(0, 5); % return convergence tracking information
    itr        = 0;           % iteration counter
    lambda     = 0;           % linear search step size
    max_affine = 60;          % maximum number of iterations in affine sub-problem
    alpha      = 0;           % initial exact penalty function parameter


    %
    % initialize current value
    xk1      = zero_n;
    %dist_cur = 0;
    for k  = 1 : N
        xk                           = x_cur(:, k);
        [gk, Gk]                     = feval(g_fun, k, xk1);
        [hk, Hk]                     = feval(h_fun, k, xk);
        %
        g_cur(:, k)                  = gk - xk;
        h_cur(:, k)                  = hk;
        %
        dg_cur(:,:, k)               = Gk;
        dh_cur(:,:, k)               = Hk;
        %
        xk1                          = xk;
    end


    % algorithm iteration loop
    for itr = 0 : max_itr

        % objective funciton corresponding to x_cur
        obj_cur = ckbs_L2L1_obj(zero_n_N, z, g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);

        % solve affine subproblem
        [x_new, r_new, s_new, p_plus_new, p_minus_new, info_affine] = ckbs_L1_affine(max_affine, delta, ...
            z, g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);

        % take the last row information from the info_affine struct and
        % append it
        info_itr = info_affine(end, :);
        info = [info; info_itr];

        if( any( any( r_new < 0 ) ) ) ...
                error('L1_nonlinear:  r_new has a negative element');
        end
        if( any( any( s_new < 0 ) ) ) ...
                error('L1_nonlinear:  s_new has a negative element');
        end
        if( any( any( p_plus_new < 0 ) ) ) ...
                error('L1_nonlinear:  p_plus_new has a negative element');
        end
        if( any( any( p_minus_new < 0 ) ) ) ...
                error('L1_nonlinear:  p_minus_new has a negative element');
        end


        % objective corresponding to x_new = x_cur + dx_cur
        obj_aff = ckbs_L2L1_obj(x_new, z, g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);


        % Compute h(F'(x)d + F(x)) - h(F(x)).
        dirDer = obj_aff - obj_cur;

        converge = (abs(dirDer) < 1e-10);
        %
        if (itr == max_itr) || (converge)
            x_out   = x_cur;
            r_out = r_new;
            s_out = s_new;
            p_plus_out = p_plus_new;
            p_minus_out = p_minus_new;

            return
        end
        if(dirDer >= 0)

            error('dirDer >= 0 in L1_nonlinear');
        end

        % ---------------------------------------------------------------
        % Line Search
        % ---------------------------------------------------------------


        c = 0.001;
        gamma = 0.5;
        lambda = 2.;


        done = false;
        search_itr = 0;
        max_search_itr = 50;

        while(~done)&&(search_itr < max_search_itr)
            search_itr = search_itr + 1;
            lambda = lambda * gamma;

            x_lambda = x_cur + lambda * (x_new);

            xk1      = zero_n;
            for k      = 1 : N
                xk                 = x_lambda(:,k);
                %
                [gk, Gk]           = feval(g_fun, k, xk1);
                [hk, Hk]           = feval(h_fun, k, xk);
                %
                %
                g_cur(:, k)        = gk - xk;
                h_cur(:,k)         = hk;
                %
                dg_cur(:,:, k)     = Gk;
                dh_cur(:,:, k)     = Hk;
                %
                xk1                = xk;
            end

            obj_lambda = ckbs_L2L1_obj(zero_n_N, z, ...
                g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);

            done = ((obj_lambda - obj_cur) <= c*lambda*dirDer);
        end

        if(search_itr == max_search_itr)
            error('line search did not converge');
        end

        x_cur = x_lambda;



    end
    error('L1_nonlinear: program error');
    return
end
