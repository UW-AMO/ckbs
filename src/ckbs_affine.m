% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Y. Aravkin: sasha.aravkin at gmail dot com
%          James V. Burke:       burke at math dot washington dot edu
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_affine$$ $newlinech %$$
% $spell
%       itr
%       complementarity
%       obj
%       ckbs
%       dg
%       dh
%       qinv
%       rinv
% $$
%
% $index ckbs_affine$$
%
% $index affine, constrained smoother$$
% $index constrained, affine smoother$$
% $index smoother, affine constrained$$
%
% $section Constrained Affine Kalman Bucy Smoother$$
%
% $head Syntax$$
% $codei/[/xOut/, /uOut/, /info/] = ckbs_affine(/max_itr/, /epsilon/, /
%       z/, /b/, /g/, /h/, /db/, /dg/, /dh/, /qinv/, /rinv/)/$$
%
% $head Purpose$$
% This routine minimizes the
% affine Kalman-Bucy smoother residual sum of squares objective
% subject to an affine inequality constraint.
%
% $head Notation$$
% The affine Kalman-Bucy smoother residual sum of squares is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = &
% \frac{1}{2}
% ( z_k - h_k - H_k * x_k )^\R{T} * R_k^{-1} * ( z_k - h_k - H_k * x_k )
% \\
% & + &
% \frac{1}{2}
% ( x_k - g_k - G_k * x_{k-1} )^\R{T} * Q_k^{-1} * ( x_k - g_k - G_k * x_{k-1} )
% \\
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
% Note that $latex g_1$$ is the initial state estimate
% and $latex Q_1$$ is the corresponding covariance.
%
% $head Problem$$
% The affine constrained Kalman-Bucy smoother problem is
% $latex \[
% \begin{array}{rll}
% {\rm minimize}
%       & S( x_1 , \ldots , x_N )
%       & {\rm w.r.t.} \; x_1 \in \B{R}^n , \ldots , x_N \in \B{R}^n
% \\
% {\rm subject \; to}
%       & b_k + B_k * x_k  \leq 0
%       & {\rm for} \; k = 1 , \ldots , N
% \end{array}
% \] $$
%
% $head First Order Conditions$$
% A state sequence $latex ( x_1 , \ldots , x_N )$$ is considered a solution
% if there is a Lagrange multiplier sequence $latex ( u_1 , \ldots , u_N )$$
% such that the following conditions are satisfied.
% $latex \[
% \begin{array}{rcl}
%       b_k + B_k * x_k                 & \leq & \varepsilon    \\
%       0                               & \leq & u_k         \\
%       | ( B_k^T * u_k + d_k )_j |        & \leq & \varepsilon \\
%       | (u_k)_i * ( b_k + B_k * x_k)_i | & \leq & \varepsilon
% \end{array}
% \] $$
% for $latex j = 1 , \ldots , n$$,
% $latex i = 1 , \ldots , \ell$$, and
% $latex k = 1 , \ldots , N$$.
% Here
% $latex d_k$$ is the partial derivative of $latex S ( x_1 , \ldots , x_N )$$
% with respect to $latex x_k$$
% and $latex (u_k)_i$$ denotes the $th i$$ component of $latex u_k$$.
%
% $head max_itr$$
% The integer scalar $icode max_itr$$ specifies the maximum number of
% iterations of the algorithm to execute. It must be greater than or
% equal to zero. Note that if it is zero, the first row of the
% $cref/info/ckbs_affine/info/$$ return value will still be computed.
% This can be useful for deciding what is a good value for the argument
% $cref/epsilon/ckbs_affine/epsilon/$$.
%
% $head epsilon$$
% The positive scalar $icode epsilon$$ specifies the convergence
% criteria value; i.e.,
% $latex \[
%       \varepsilon = epsilon
% \] $$
%
% $head z$$
% The argument $icode z$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       z_k = z(:, k)
% \]$$
% and $icode z$$ has size $latex m \times N$$.
%
% $head b$$
% The argument $icode b$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       b_k = b(:, k)
% \]$$
% and $icode b$$ has size $latex \ell \times N$$.
% If $latex \ell = 0$$, the problem is not constrained; i.e.,
% it is the affine Kalman-Bucy smoother problem.
%
% $head g$$
% The argument $icode g$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       g_k = g(:, k)
% \]$$
% and $icode g$$ has size $latex n \times N$$.
% The value $latex g_1$$ serves as the initial state estimate.
%
% $head h$$
% The argument $icode h$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       h_k = h(:, k)
% \]$$
% and $icode h$$ has size $latex m \times N$$.
%
% $head db$$
% The argument $icode db$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       B_k = db(:,:,k)
% \]$$
% and $icode db$$ has size $latex \ell \times n \times N$$.
%
%
% $head dg$$
% The argument $icode dg$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       G_k = dg(:,:,k)
% \]$$
% and $icode dg$$ has size $latex n \times n \times N$$.
% The initial state estimate $latex g_1$$ does not depend on the value of
% $latex x_0$$, hence $latex G_1$$ should be zero.
%
% $head dh$$
% The argument $icode dh$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       H_k = dh(:,:,k)
% \]$$
% and $icode dh$$ has size $latex m \times n \times N$$.
%
% $head qinv$$
% The argument $icode qinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       Q_k^{-1} = qinv(:,:, k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
% The value of $latex Q_k$$ is the variance of the initial state
% estimate $latex g_1$$.
%
% $head rinv$$
% The argument $icode rinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       R_k^{-1} = rinv(:,:, k)
% \]$$
% and $icode rinv$$ has size $latex m \times m \times N$$.
% It is ok to signify a missing data value by setting the corresponding
% row and column of $icode rinv$$ to zero. This also enables use
% to interpolate the state vector $latex x_k$$ to points where
% there are no measurements.
%
% $head xOut$$
% The result $icode xOut$$ contains the optimal sequence
% $latex ( x_1 , \ldots , x_N )$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       x_k = xOut(:, k)
% \]$$
% and $icode xOut$$ is a two dimensional array with size $latex n \times N$$.
%
% $head uOut$$
% The result $icode uOut$$ contains the Lagrange multiplier sequence
% corresponding to $icode xOut$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       u_k = uOut(:, k)
% \]$$
% and $icode uOut$$ is a two dimensional array with size
% $latex \ell \times N$$.
% The pair $icode xOut$$, $icode uOut$$ satisfy the
% $cref/first order conditions/ckbs_affine/First Order Conditions/$$.
%
% $head info$$
% The result $icode info$$ is a matrix with each row corresponding
% to an iteration of the algorithm.
% Note that $code ckbs_affine$$ has satisfied the convergence condition if
% and only if
% $codei%
%       all( %info%(end, 1:3) <= %epsilon% )
% %$$
% $pre
%
% $$
% We use $latex x_k^q$$ ($latex u_k^q$$)
% to denote the state vector (dual vector)
% for time point $latex k$$ and at the end of iteration $latex q-1$$.
% We use $latex d_k^q$$ for the partial derivative of
% $latex S ( x_1^q , \ldots , x_N^q )$$ with respect to $latex x_k$$.
%
% $subhead Constraint Bound$$
% The value $icode%info%(%q%, 1)%$$ is
% a bound on the constraint functions
% at the end of iteration $latex q-1$$. To be specific
% $latex \[
% info(q, 1) = \max_{i,k} ( b_k + B_k * x_k )_i
% \] $$
%
% $subhead First Order Conditions$$
% The value $icode%info%(%q%, 2)%$$ is
% a bound on the partial derivative of the Lagrangian with respect to
% the state vector sequence
% at the end of iteration $latex q-1$$:
% $latex \[
% info(q, 2) = \max \left[ | ( B_k^T * u_k + d_k )_j | \right]
% \] $$
% where the maximum is with respect to $latex j = 1 , \ldots , n$$
% and $latex k = 1 , \ldots , N$$.
%
% $subhead Complementarity Conditions$$
% The value $icode%info%(%q%, 3)%$$ is
% a bound on the complementarity conditions
% at the end of iteration $latex q-1$$:
% $latex \[
% info(q, 3) = \max \left[ | (u_k)_i * ( b_k + B_k * x_k)_i | \right]
% \] $$
% where the maximum is with respect to
% $latex i = 1 , \ldots , \ell$$ and
% $latex k = 1 , \ldots , N$$.
%
%
% $subhead Step Size$$
% The value $icode%info%(%q%, 4)%$$ is the line search step size used during
% iteration $latex q-1$$.
% Small step sizes indicate problems with the interior point method
% used to solve the affine problem
% (with the exception that $icode%info%(1, 4)%$$ is always zero).
%
%
% $children#
%       example/affine_ok_box.m
% #$$
%
% $head Example$$
% The file $cref affine_ok_box.m$$ contains an example and test of
% $code ckbs_affine$$.
% It returns true if $code ckbs_affine$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------
function [xOut, uOut, info] = ckbs_affine(max_itr, epsilon, ...
    z, b, g, h, db, dg, dh, qinv, rinv)
    if nargin ~= 11
        error('ckbs_affine: improper number of input arguments');
    end
    if nargout ~= 3
        error('ckbs_affine: improper number of return values');
    end
    % size of problem
    n     = size(g, 1);
    N     = size(z, 2);
    m     = size(z,   1);
    ell   = size(b,   1);
    %
    % check sizes
    if N ~= size(b,2) | N ~= size(h,2) | N ~= size(db,3) | ...
            N ~= size(dg,3) | N ~= size(dh,3) | N ~= size(qinv,3) | N ~= size(rinv,3)
        N
        size(z,2)
        size(b,2)
        size(h,2)
        size(db,3)
        size(dg,3)
        size(dh,3)
        size(qinv,3)
        size(rinv,3)
        error('ckbs_affine: argument sizes do not agree');
    end
    if n ~= size(db,2) | n ~= size(dg,2) | n ~= size(dh,2) | ...
            n ~= size(qinv,1) | n ~= size(qinv,2)
        n
        size(g,1)
        size(db,2)
        size(dg,2)
        size(dh,2)
        size(qinv,1)
        size(qinv,2)
        error('ckbs_affine: argument sizes do not agree');
    end
    if m ~= size(h,1) | m ~= size(dh,1) | m ~= size(rinv,1) | m ~= size(rinv,2)
        m
        size(h,1)
        size(dh,1)
        size(rinv,1)
        size(rinv,2)
        error('ckbs_affine: argument sizes do not agree');
    end
    if ell ~= size(db,1)
        ell
        size(db,1)
        error('ckbs_affine: argument sizes do not agree');
    end
    %
    % Other usefull sizes
    r        = ell * N;
    p        = n   * N;
    %
    % vector of zero state values
    xZero = zeros(n, N);
    %
    % diagonal and off diagonal blocks for Hessian of the objective
    [D, A] = ckbs_sumsq_hes(dg, dh, qinv, rinv);
    %
    % gradient of the objective that corresponds to xIn
    d      = ckbs_sumsq_grad(xZero, z, g, h, dg, dh, qinv, rinv);
    dVec   = reshape(d, p, 1);
    %
    % solve the unconstrained problem
    
%     fid = fopen('BadDiag.dat', 'wt');
%     fprintf(fid, '%6.6f %6.6f\n', D);
%     fid = fopen('BadOffDiag.dat', 'wt');
%     fprintf(fid, '%6.6f %6.6f\n', A);
%     fid = fopen('BadVec.dat', 'wt');
%     fprintf(fid, '%6.6f %6.6f\n', -dVec);
    
    [y, lambda] = ckbs_tridiag_solve(D, A, - dVec);
    xOut        = reshape(y, n, N);
    uOut        = zeros(ell, N);
    %
    % gradient at solution of unconstrained problem
    d          = ckbs_sumsq_grad(xOut, z, g, h, dg, dh, qinv, rinv);
    e2         = max( max( abs(d) ) );
    if e2 > epsilon
        'warning: ckbs_affine: cant solve unconstrained problem within epsilon'
        e2
        epsilon
    end
    %
    % check for no constraints
    if ell == 0
        % check if gradient is as small as requested
        info.iters     = [ 0, e2 , 0 , 0 ];
        info.D = D;
        info.A = A;
        return
    end
    %
    % vector version of b
    bVec   = reshape(b, r, 1);
    %
    % check if unconstrained solution solve convergence criteria
    bTmp = bVec + ckbs_blkdiag_mul(db, y);
    e1   = max( bTmp );
    info.iters = [ e1, e2 , 0 , 0 ];
    if e1 <= epsilon
        info.D = D;
        info.A = A;
        return;
    end
    %
    % initialiize mu, s, u to be near central path corresponding to y
    mu = e1;
    s  = sqrt(mu) * ones(r, 1);
    u  = sqrt(mu) * ones(r, 1);
    %
    if min(s) <= 0
        error('ckbs_affine: initail s not all positive');
    end
    %
    % how close to boundry are we willing to go in one step
    gamma   = .1;
    %
    % determine the value of y that solves the problem
    converge = false;
    itr       = 0;
    F         = ckbs_kuhn_tucker(mu, s, y, u, bVec, dVec, db, D, A);
    while ( ~ converge ) & (itr < max_itr)
        itr = itr + 1;
        %
        % Newton step corresponding to the mu-penalaty problem
        [ds, dy, du] = ckbs_newton_step(mu, s, y, u, bVec, dVec, db, D, A);
        %
        % determine maximum allowable step factor lambda
        ratio     = [ ds ; du ] ./ [ s ; u ];
        ratio_max = max( - ratio );
        if ratio_max >= 1 - gamma
            lambda = (1 - gamma) / ratio_max;
        else
            lambda = 1;
        end
        %
        % line search
        %
        ok        = 0;
        kount     = 0;
        max_kount = 35;
        lambda    = 2 * lambda;
        while (~ok) & (kount < max_kount)
            kount  = kount + 1;
            lambda = lambda / 2;
            %
            % step of size lambda
            s_new = s + lambda * ds;
            y_new = y + lambda * dy;
            u_new = u + lambda * du;
            %
            % check for feasibility
            if min(s_new) <= 0 | min(u_new) <= 0
                error('ckbs_affine: program error');
            end
            F_new = ckbs_kuhn_tucker( ...
                mu, s_new, y_new, u_new, bVec, dVec, db, D, A ...
                );
            % Can use any norm here
            % G     = sqrt( F' * F );
            % G_new = sqrt( F_new' * F_new );
            G     = max( abs( F ) );
            G_new = max( abs( F_new ) );
            ok    = G_new <= (1 - gamma * lambda ) * G;
        end
        lambda = lambda;
        if ~ok
            error('ckbs_affine: linear search did not converge');
        end
        F    = F_new;
        s    = s_new;
        y    = y_new;
        u    = u_new;
        %
        % new value for objective function
        xOut     = reshape(y, n, N);
        uOut     = reshape(u, ell, N);
        %
        % check for convergence
        %Bt_u_d   = F((r+1) : (r+p));
        e1       = norm(F(1:r), inf);
        e2       = norm(F(r+1:r+p), inf);
        e3       = norm(u .* s, inf );
        info_itr = [ e1, e2, e3, lambda ];
        info.iters     = [ info.iters ; info_itr ];
        converge = (e1 < epsilon) & (e2 < epsilon) & (e3 <= epsilon);
        if converge || kount == max_kount
            info.D = D;
            info.A = A;
        end
        %
        % every third step is a corrector
        if ( mod(itr, 3) ~= 1 )
            %mu_old = mu;
            %mu     = mu / 10;

            comp = max(u.*s);
            compFrac = comp/(m*N);
            mu = compFrac*0.1;
            %dmu    = mu - mu_old;
            %
            % Kuhn-Tucker residuals for new relaxation parameter
            F = ckbs_kuhn_tucker( ...
                mu, s, y, u, bVec, dVec, db, D, A ...
                );
            %F( (r+p+1) : (r+p+r) ) = F( (r+p+1) : (r+p+r) ) - dmu;
        end
    end
    return
end
