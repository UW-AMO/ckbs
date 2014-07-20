% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_sumsq_grad$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_sumsq_grad
%       dg
%       dh
%       qinv
%       rinv
% $$
%
% $section Affine Residual Sum of Squares Gradient$$
%
% $index ckbs_sumsq_grad$$
% $index sumsq_grad$$
%
% $index gradient, of objective$$
% $index objective, gradient$$
% $index residual, gradient$$
% $index sum of squares, gradient$$
%
% $head Syntax$$
% $codei/[/grad/] = ckbs_sumsq_grad(/
%       x/, /z/, /g/, /h/, /dg/, /dh/, /qinv/, /rinv/)/$$
%
% $head Purpose$$
% This computes the gradient of the
% of the affine Kalman-Bucy smoother residual sum of squares.
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
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
%
% $head Gradient$$
% We define $latex Q_{N+1}$$ to be the $latex n \times n$$ identity
% matrix and $latex G_{N+1}$$ to be zero,
% $latex \[
% \begin{array}{rcl}
% \nabla_k S_k^{(1)} ( x_k , x_{k-1} )
% & = &  H_k^\R{T} * R_k^{-1} * ( h_k + H_k * x_k - z_k )
%   +    Q_k^{-1} * ( x_k - g_k - G_k * x_{k-1} )
% \\
% \nabla_k S_{k+1}^{(1)} ( x_{k+1} , x_k )
% & = & G_{k+1}^\R{T} * Q_{k+1}^{-1} * ( g_{k+1} + G_{k+1} * x_k  - x_{k+1} )
% \end{array}
% \] $$
% It follows that the gradient of the
% affine Kalman-Bucy smoother residual sum of squares is
% $latex \[
% \begin{array}{rcl}
% \nabla S ( x_1 , \ldots , x_N )
% & = &
% \left( \begin{array}{c}
%       d_1 \\ \vdots \\ d_N
% \end{array} \right)
% \\
% d_k & = & \nabla_k S_k^{(1)}     ( x_k , x_{k-1} )
%       +   \nabla_k S_{k+1}^{(1)} ( x_{k+1} , x_k )
% \end{array}
% \] $$
% where $latex S_{N+1} ( x_{N+1} , x_N )$$ is defined as
% identically zero.
%
% $head x$$
% The argument $icode x$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x_k = x(:, k)
% \]$$
% and $icode x$$ has size $latex n \times N$$.
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
%
% $head h$$
% The argument $icode h$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       h_k = h(:, k)
% \]$$
% and $icode h$$ has size $latex m \times N$$.
%
% $head dg$$
% The argument $icode dg$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       G_k = dg(:,:,k)
% \]$$
% and $icode dg$$ has size $latex n \times n \times N$$.
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
%       Q_k^{-1} = qinv(:,:,k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
%
% $head rinv$$
% The argument $icode rinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       R_k^{-1} = rinv(:,:,k)
% \]$$
% and $icode rinv$$ has size $latex m \times m \times N$$.
%
% $head grad$$
% The result $icode grad$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       d_k = grad(:, k)
% \]$$
% and $icode grad$$ has size $latex n \times N$$.
%
% $children#
%       example/sumsq_grad_ok.m
% #$$
%
% $head Example$$
% The file $cref sumsq_grad_ok.m$$ contains an example and test of
% $code ckbs_sumsq_grad$$.
% It returns true if $code ckbs_sumsq_grad$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------
function grad = ckbs_sumsq_grad(x, z, g, h, dg, dh, qinv, rinv)
    %
    % sizes for this problem
    m = size(dh, 1);
    n = size(dh, 2);
    N = size(dh, 3);
    %
    % initialize return values
    grad = zeros(n, N);
    %
    xk    = zeros(n, 1);
    gk    = zeros(n, 1);
    Gk    = zeros(n, n);
    Qkinv = eye(n, n);
    for k = N : -1 : 1
        if k > 1
            xkm = x(:, (k-1));
        else
            xkm = zeros(n, 1);
        end
        xk1     = xk;
        gk1     = gk;
        Gk1     = Gk;
        Qk1inv  = Qkinv;
        %
        xk      = x(:, k);
        zk      = z(:, k);
        hk      = h(:, k);
        gk      = g(:, k);
        Gk      = dg(:,:, k);
        Qkinv   = qinv(:,:, k);
        Hk      = dh(:,:, k);
        Rkinv   = rinv(:,:, k);
        %
        xk1res  = gk1 + Gk1 * xk - xk1;
        xkres   = xk - gk - Gk * xkm;
        zkres   = hk + Hk * xk - zk;
        %
        grad(:, k) = Hk' * Rkinv * zkres ...
            + Qkinv * xkres ...
            + Gk1'* Qk1inv * xk1res;
    end
    return
end
