% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_L1_affine$$ $newlinech %$$
% $spell
%       itr
%       complementarity
%       obj
%       ckbs
%       dg
%       dh
%       qinv
%       rinv
%       optimality
%       pairwise
%       Kuhn
% $$
%
% $index ckbs_L1_affine$$
%
% $index affine, robust smoother$$
% $index robust, affine smoother$$
% $index smoother, affine robust$$
%
% $section Robust L1 Affine Kalman Bucy Smoother$$
%
% $head Syntax$$
% $codei/[/xOut/, /rOut/, /sOut/, /pPlusOut/, /pMinusOut/, /info/] = ckbs_L1_affine(...
%        /max_itr/, /epsilon/, /z/, /b/, /g/, /h/, /db/, ...
%        /dg/, /dh/, /qinv/, /rinv/)/$$
%
% $head Purpose$$
% This routine minimizes the
% robust L1 Kalman-Bucy smoother objective for affine process and
% measurement models.
%
% $head Notation$$
% The robust L1 Kalman-Bucy smoother objective is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = &
% \frac{1}{2}
% \|R_k^{-1/2}( z_k - h_k - H_k * x_k )^\R{T}\|_1
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
% The robust L1 Kalman-Bucy smoother problem is
%$latex \[
%\begin{array}{ll}
%{\rm minimize}
%& \frac{1}{2} y^\R{T} H(x) y + d(x)^\R{T} y
%       + \sqrt{\B{2}}^\R{T} (p^+ + p^-)
%    \;{\rm w.r.t} \;  y \in \B{R}^{nN}\; , \; p^+ \in \B{R}_+^{M} \; , \; p^- \in \B{R}_+^{M}
%\\
%{\rm subject \; to} &  b(x) + B(x) y - p^+ + p^- = 0
%\end{array}
%\] $$
%
% $head Lagrangian and First Order Conditions$$
% We denote the lagrange multiplier associated to the
% equality constraint by $latex q$$, the lagrange multiplier
% associated to the non-negativity constraint on $latex p^+$$ by
% $latex r$$, and the lagrange multiplier associated to the
% non-negativity of $latex p^-$$ by $latex s$$. We also use  $latex \B{1}$$
% to denote the vector of length $latex mN$$ with all its components
% equal to one, and $latex \B{\sqrt{2}}$$ to denote the vector of
% length $latex mN$$ with all its components equal to $latex
% \sqrt{2}$$.
% The Lagrangian corresponding to the L1 Kalman-Bucy problem is
% given by
% $latex \[
% L(p^+, p^-, y, q, r, s)  =
% \frac{1}{2} y^\R{T} H y + d^\R{T} y + \B{\sqrt{2}}^T(p^+ + p^-)
% + q^\R{T} (b + B y - p^+ + p^-) - r^\R{T} p^+ - s^\R{T} p^- \;.
% \] $$
% The partial gradients of the Lagrangian are given by
% $latex \[
% \nabla L :=
% \begin{array}{rcl}
% \nabla_p^+ L(p^+, p^-, y, q ) & = & \B{\sqrt{2}} - q - r \\
% \nabla_p^- L(p^+, p^-, y, q) & = & \B{\sqrt{2}} + q - s \\
% \nabla_y L(p^+, p^-, y, q ) & = & H y + c + B^\R{T} q \\
% \nabla_q L(p^+, p^-, y, q ) & = & b + B y - p^+ + p^- \\
% \end{array}
% \] $$
% At optimality, all of these gradients will be zero
% (or smaller than a tolerance $latex \varepsilon$$), and
% complementarity conditions hold, i.e.
% $latex
% \[ r_i p_i^+ = s_i p_i^- = 0 \; \forall \; i \;.
% \]
% $$ (or again, each pairwise product is smaller than $latex
% \varepsilon$$.
%
%
% $head max_itr$$
% The integer scalar $icode max_itr$$ specifies the maximum number of
% iterations of the algorithm to execute. It must be greater than or
% equal to zero.
%
% $head epsilon$$
% The positive scalar $icode epsilon$$ specifies the convergence
% criteria value; i.e.,
% $latex \[
%       \|\nabla L\|_1 \leq epsilon\;.
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
% $head rOut$$
% The result $icode rOut$$ contains the Lagrange multiplier sequence
% corresponding to the non-negativity constraint on $latex p^+$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       r_k = rOut(:, k)
% \]$$
% and $icode rOut$$ is a two dimensional array with size
% $latex m \times N$$.
%
% $head sOut$$
% The result $icode rOut$$ contains the Lagrange multiplier sequence
% corresponding to the non-negativity constraint on $latex p^-$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       s_k = sOut(:, k)
% \]$$
% and $icode sOut$$ is a two dimensional array with size
% $latex m \times N$$.
%
% $head pPlusOut$$
% The result $icode pPlusOut$$ contains the positive parts
% of the measurement residual $latex b(x) + B(x) y$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       p_k = pPlusOut(:, k)
% \]$$
% and $icode pPlusOut$$ is a two dimensional array with size
% $latex m \times N$$.
%
% $head pMinusOut$$
% The result $icode pMinusOut$$ contains the negative parts
% of the measurement residual $latex b(x) + B(x) y$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       p_k = pMinusOut(:, k)
% \]$$
% and $icode pMinusOut$$ is a two dimensional array with size
% $latex m \times N$$.
%
% $head info$$
% The result $icode info$$ is a matrix with each row corresponding
% to an iteration of the algorithm.
% There are five numbers output per iteration:
% the infinity norm of the Kuhn-Tucker conditions,
% the one-norm of the Kuhn-Tucker conditions, the duality gap,
% $latex \mu$$, and the number of line search
% Note that $code ckbs_L1_affine$$ has satisfied the convergence condition if
% $codei%
%        %info%(end, 2) <= %epsilon%
% %$$
% $pre
%
% $$
%
% $children#
%       example/L1_affine_ok.m
% #$$
%
% $head Example$$
% The file $cref L1_affine_ok.m$$ contains an example and test of
% $code ckbs_L1_affine$$.
% It returns true if $code ckbs_L1_affine$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------


function [xOut, rOut, sOut, pPlusOut, pMinusOut, info] = ckbs_L1_affine(max_itr, epsilon, z, g, h, dg, dh, qinv, rinv)

    info = [];
    if nargin ~= 9
        error('ckbs_affine: improper number of input arguments');
    end
    if nargout ~= 6
        error('ckbs_affine: improper number of return values');
    end
    % size of problem
    n     = size(g, 1);
    N     = size(z, 2);
    m     = size(z,   1);
    %ell   = size(b,   1);
    %
    % check sizes
    if  N ~= size(h,2) || N ~= size(dg,3) || N ~= size(dh,3) || N ~= size(qinv,3) || N ~= size(rinv,3)
        N
        size(z,2)
        size(h,2)
        size(dg,3)
        size(dh,3)
        size(qinv,3)
        size(rinv,3)
        error('ckbs_affine: argument sizes do not agree');
    end
    if n ~= size(dg,2) || n ~= size(dh,2) || ...
            n ~= size(qinv,1) ||n ~= size(qinv,2)
        n
        size(g,1)
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
    %
    % Other usefull sizes
    %meas        = ell * N;
    obs        = n   * N;
    %
    % vector of zero state values
    xZero = zeros(n, N);


    % Update D and A to not include measurements
    [D, A] = ckbs_process_hes(dg, qinv);

    a = ckbs_process_grad(xZero, g, dg, qinv);
    c = a(:);



    y = ckbs_tridiag_solve_b(D, A, -c(:));

    cRec = ckbs_blktridiag_mul(D, A, y);

    if(norm(-c(:) - cRec(:), inf) > 1e-5)
        norm(-c(:) - cRec(:))
        error('tridiagonal solver failed');
    end



    y = reshape(y, n, N);


    diff = zeros(m, N);
    b = zeros(m, N);
    r = zeros(m, N);
    s = zeros(m, N);


    crinv = zeros(m, m, N);     % Cholesky decomposition of rinv
    B = zeros(m, n, N);
    pPlus = zeros(m, N);        % max(0, b + By)
    pMinus = zeros(m,N);        % max(0, -b + BY)



    for k = 1:N


        crinv(:,:,k) = sqrt(rinv(:,:,k));
        B(:,:,k) = -crinv(:,:,k)*dh(:,:,k);
        b(:,k) = crinv(:,:,k)*(z(:,k) - h(:,k));

        %Version 1
        temp = b(:,k) + B(:,:,k)*y(:,k);
        pPlus(:,k) = 500 + max(0, temp); % initialization:
        pMinus(:,k) = 500 + max(0, -temp);
    end

    mu = 100;

    for k = 1:N

        r(:,k) = sqrt(2)*ones(m,1);
        s(:,k) = sqrt(2)*ones(m,1);

    end

    if (min(min(s)) <= 0) || (min(min(r)) <=0) || (min(min(pPlus)) <=0) || (min(min(pMinus)) <=0)
        error('L1_affine: initial s, r, pPlus, pMinus not all positive');
    end

    yVec = y; % in case we need this
    y = reshape(yVec, n, N); % I assume y is in this form

    % how close to boundry are we willing to go in one step
    gamma   = .01;
    %
    % determine the value of y that solves the problem
    converge = false;
    itr       = 0;



    while ( ~ converge ) && (itr < max_itr)



        itr = itr + 1;
        F         = ckbs_kuhn_tucker_L1(mu, s, y, r, b, a, B, D, A,  pPlus, pMinus);

        [dPPlus, dPMinus, dr, ds, dy] = ckbs_newton_step_L1(mu, s, y, r, b, a, B, D, A, pPlus, pMinus);



        % SHOULD BE 0!!!
        expr1 =    norm(pMinus(:).*s(:) + s(:).*dPMinus(:) + pMinus(:).*ds(:) - mu, inf);
        expr2 =    norm(pPlus(:).*r(:) + r(:).*dPPlus(:) + pPlus(:).*dr(:) - mu, inf);

        if (abs(expr1) > 1e-4) || (abs(expr2) > 1e-4)

            expr1
            expr2
            error('L1`_affine: Newton Solver not working, expr not 0');


        end


        %
        % determine maximum allowable step factor lambda
        ratio     = [ dr ; ds; dPPlus; dPMinus ] ./ [ r ; s; pPlus; pMinus ];



        ratioMax = max(max( - ratio ));

        if (ratioMax <=0)
            lambda = 1;

        else
            rNeg = -1./ratio(ratio < 0);
            %min(min(ratio))
            maxNeg = min(min(rNeg));
            lambda = .99*min(maxNeg, 1);
        end


        % line search
        %
        ok        = 0;
        kount     = 0;
        max_kount = 18;
        beta = 0.5;
        lambda = lambda/beta;
        while (~ok) && (kount < max_kount)
            kount  = kount + 1;


            lambda = lambda *beta;
            %
            % step of size lambda
            s_new = s + lambda * ds;
            y_new = y + lambda * dy;
            r_new = r + lambda * dr;
            pPlus_new = pPlus + lambda*dPPlus;
            pMinus_new = pMinus + lambda*dPMinus;
            %


            % check for feasibility
            if min(min(s_new)) <= 0 || min(min(r_new)) <=0 ||min(min(pPlus_new))<=0 || min(min(pMinus_new)) <= 0
                error('L1_affine: program error, negative entries');
            end


            F_new = ckbs_kuhn_tucker_L1(...
                mu, s_new, y_new, r_new, b, a, B, D, A, pPlus_new, pMinus_new...
                );


            G     = max(abs(F));
            G_new = max(abs(F_new));

            ok   = (G_new <= (1 - gamma *lambda) * G);
        end

        if ~ok
            df = max(F - F_new);
            if(df <= epsilon)
                return
            end
            warning('L1_affine: line search failed');
        end

        F    = F_new;
        %
        s    = s_new;
        y    = y_new;
        r    = r_new;
        pPlus = pPlus_new;
        pMinus = pMinus_new;
        %
        % new value for objective function
        xOut     = y;
        rOut     = r;
        sOut     = s;
        pPlusOut = pPlus;
        pMinusOut = pMinus;
        muOut = mu;

        VP = ckbs_L2L1_obj(y, z, g, h, dg, dh, qinv, rinv);
        Kxnu = ckbs_L2L1_obj(xZero, z, g, h, dg, dh, qinv, rinv);

        G1 = sum(sum(r.*pPlus)) + sum(sum(s.*pMinus));
        converge = (G1 < min(Kxnu - VP, epsilon));
        %converge = (G1 < Kxnu - VP);

        % every third step is a corrector
        if ( mod(itr, 3) ~= 1 )

            temp = sum(r(:).*pPlus(:)) + sum(s(:).*pMinus(:));
            compMuFrac = temp/(2*m*N);
            muNew = .1*compMuFrac;
            mu = muNew;


        end
        info = [info; G, G1, Kxnu - VP, mu, kount];
    end
    return
end
