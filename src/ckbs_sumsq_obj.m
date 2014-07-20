% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at us dot ibm dot com
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_sumsq_obj$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_sumsq_obj
%       qinv
%       rinv
%       dh
%       dg
% $$
%
% $section Affine Residual Sum of Squares Objective$$
%
% $index ckbs_sumsq_obj$$
% $index sumsq_obj$$
%
% $index objective, residual sum of squares$$
% $index residual, sum of squares objective$$
% $index sum of squares, residual objective$$
%
% $head Syntax$$
% $codei/[/obj/] = ckbs_sumsq_obj(/
%       x/, /z/, /g/, /h/, /dg/, /dh/, /qinv/, /rinv/)/$$
%
% $head Purpose$$
% This routine computes the value of the
% affine Kalman-Bucy smoother residual sum of squares objective function.
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
%       Q_k^{-1} = qinv(:,:, k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
%
% $head rinv$$
% The argument $icode rinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       R_k^{-1} = rinv(:,:, k)
% \]$$
% and $icode rinv$$ has size $latex m \times m \times N$$.
%
% $head obj$$
% The result $icode obj$$ is a scalar given by
% $latex \[
%       obj = S ( x_1 , \ldots , x_N )
% \] $$
%
% $children#
%       example/sumsq_obj_ok.m
% #$$
%
% $head Example$$
% The file $cref sumsq_obj_ok.m$$ contains an example and test of
% $code ckbs_sumsq_obj$$.
% It returns true if $code ckbs_sumsq_obj$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------
function [obj] = ckbs_sumsq_obj(x, z, g, h, dg, dh, qinv, rinv)
    n    = size(x, 1);
    N    = size(x, 2);
    obj  = 0;
    xk   = zeros(n, 1);
    for k = 1 : N
        xkm   = xk;
        xk    = x(:, k);
        zk    = z(:, k);
        hk    = h(:, k);
        gk    = g(:, k);
        Hk    = dh(:,:, k);
        Gk    = dg(:,:, k);
        Qkinv = qinv(:, :, k);
        Rkinv = rinv(:, :, k);
        zres  = zk - hk - Hk * xk;
        xres  = xk - gk - Gk * xkm;
        sk    = (zres' * Rkinv * zres + xres' * Qkinv * xres) / 2;
        obj   = obj + sk;
    end
    return
end
