% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_sumsq_hes$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_sumsq_hes
%       dg
%       dh
%       tri
%       qinv
%       rinv
% $$
%
% $section Affine Residual Sum of Squares Hessian$$
%
% $index ckbs_sumsq_hes$$
% $index sumsq_hes$$
%
% $index Hessian, of objective$$
% $index objective, Hessian$$
% $index residual, Hessian$$
% $index sum of squares, Hessian$$
%
% $head Syntax$$
% $codei/[/D/, /A/] = ckbs_sumsq_hes(/dg/, /dh/, /qinv/, /rinv/)/$$
%
% $head Purpose$$
% This routine returns the diagonal and off diagonal blocks corresponding to
% the Hessian of the affine Kalman-Bucy smoother residual sum of squares.
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
% $head Hessian$$
% If we define $latex Q_{N+1}$$ to be the $latex n \times n$$ identity
% matrix and $latex G_{N+1}$$ to be zero,
% the Hessian of the
% affine Kalman-Bucy smoother residual sum of squares is
% $latex \[
% \begin{array}{rcl}
% S^{(2)} ( x_1 , \ldots , x_N ) & = &
% \left( \begin{array}{cccc}
% D_1 & A_2^\R{T} & 0         &           \\
% A_2 & D_2       & A_3^\R{T} & 0         \\
% 0   & \ddots    & \ddots    & \ddots    \\
%     & 0         & A_N       & D_N
% \end{array} \right)
% \\
% D_k & = & H_k^\R{T} * R_k^{-1} * H_k + Q_k^{-1}
%       + G_{k+1}^\R{T} * Q_{k+1}^{-1} * G_{k+1}
% \\
% A_k & = & - Q_k^{-1} * G_k
% \end{array}
% \] $$
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
% $head D$$
% The result $icode D$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       D_k = D(:,:,k)
% \]$$
% and $icode D$$ has size $latex n \times n \times N$$.
%
% $head A$$
% The result $icode A$$ is a three dimensional array,
% for $latex k = 2 , \ldots , N$$
% $latex \[
%       A_k = A(:,:,k)
% \]$$
% and $icode A$$ has size $latex n \times n \times N$$.
%
%
% $children#
%       example/sumsq_hes_ok.m
% #$$
%
% $head Example$$
% The file $cref sumsq_hes_ok.m$$ contains an example and test of
% $code ckbs_sumsq_hes$$.
% It returns true if $code ckbs_sumsq_hes$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------
function [D, A] = ckbs_sumsq_hes(dg, dh, qinv, rinv)
    %
    % sizes for this problem
    m = size(dh, 1);
    n = size(dh, 2);
    N = size(dh, 3);
    %
    % check sizes
    if N ~= size(dg,3) | N ~= size(qinv,3) | N ~= size(rinv,3)
        N
        size(dg,3)
        size(qinv,3)
        size(rinv,3)
        error('ckbs_sumsq_hes: argument sizes do not agree');
    end
    if n ~= size(dg,1) | n ~= size(dg,2) | n ~= size(qinv,1) | n ~= size(qinv,2)
        n
        size(dg,1)
        size(dg,2)
        size(qinv,1)
        size(qinv,2)
        error('ckbs_sumsq_hes: argument sizes do not agree');
    end
    if m ~= size(rinv,1) | m ~= size(rinv,2)
        m
        size(rinv,1)
        size(rinv,2)
        error('ckbs_sumsq_hes: argument sizes do not agree');
    end
    %
    % dimension return values
    D  = zeros(n, n, N);
    A  = zeros(n, n, N);
    %
    % initilaize for loop
    Gk    = zeros(n, n);
    Qkinv = eye(n, n);
    for k = N : -1 : 1
        Gk1     = Gk;
        Gk      = dg(:, :, k);
        %
        Qk1inv   = Qkinv;
        Qkinv    = qinv(:, :, k);
        %
        Hk       = dh(:, :, k);
        Rkinv    = rinv(:, :, k);
        %
        D(:,:,k) = Hk' * Rkinv * Hk + Qkinv + Gk1' * Qk1inv * Gk1;
        A(:,:,k) = - Qkinv * Gk;
    end
    return
end
