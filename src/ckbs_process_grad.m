% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_process_grad$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_process_grad
%       dg
%       dh
%       qinv
%       rinv
% $$
%
% $section Affine Residual Process Sum of Squares Gradient$$
%
% $index ckbs_process_grad$$
% $index process_grad$$
%
% $index gradient, of process objective$$
% $index objective, process gradient$$
% $index process residual, gradient$$
% $index sum of squares of process, gradient$$
%
% $head Syntax$$
% $codei/[/grad/] = ckbs_process_grad(/
%       x/, /g/, /dg/, /qinv/)/$$
%
% $head Purpose$$
% This computes the gradient of the
% of the affine Kalman-Bucy smoother process residual sum of squares.
%
% $head Notation$$
% The affine Kalman-Bucy smoother process residual sum of squares is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = &
% \frac{1}{2}
% ( x_k - g_k - G_k * x_{k-1} )^\R{T} * Q_k^{-1} * ( x_k - g_k - G_k * x_{k-1} )
% \end{array}
% \] $$
% where the matrix $latex Q_k$$ is
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
%
% $head Gradient$$
% We define $latex Q_{N+1}$$ to be the $latex n \times n$$ identity
% matrix and $latex G_{N+1}$$ to be zero,
% $latex \[
% \begin{array}{rcl}
% \nabla_k S_k^{(1)} ( x_k , x_{k-1} )
% & = &  Q_k^{-1} * ( x_k - g_k - G_k * x_{k-1} )
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
% $head g$$
% The argument $icode g$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       g_k = g(:, k)
% \]$$
% and $icode g$$ has size $latex n \times N$$.
%
% $head dg$$
% The argument $icode dg$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       G_k = dg(:,:,k)
% \]$$
% and $icode dg$$ has size $latex n \times n \times N$$.
%
% $head qinv$$
% The argument $icode qinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       Q_k^{-1} = qinv(:,:,k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
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
%       example/process_grad_ok.m
% #$$
%
% $head Example$$
% The file $cref process_grad_ok.m$$ contains an example and test of
% $code ckbs_process_grad$$.
% It returns true if $code ckbs_process_grad$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------
function grad = ckbs_process_grad(x, g, dg, qinv)
    %
    % sizes for this problem
    n = size(dg, 1);
    N = size(dg, 3);
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
        gk      = g(:, k);
        Gk      = dg(:,:, k);
        Qkinv   = qinv(:,:, k);
        %
        xk1res  = gk1 + Gk1 * xk - xk1;
        xkres   = xk - gk - Gk * xkm;
        %
        grad(:, k) = Qkinv * xkres ...
            + Gk1'* Qk1inv * xk1res;
    end
    return
end
