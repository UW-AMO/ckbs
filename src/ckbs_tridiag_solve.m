% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2014
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    sasha.aravkin at gmail dot com
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_tridiag_solve$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_tridiag_solve
%       Tridiagonal
%       Bradley
% $$
%
% $section Symmetric Block Tridiagonal Algorithm$$
%
% $index ckbs_tridiag_solve$$
% $index tridiag_solve$$
%
% $index solve, block tridiagonal$$
% $index tridiagonal, block solve$$
%
% $head Syntax$$
% $codei/[/e/, /lambda/] = ckbs_tridiag_solve(/b/, /c/, /r/)/$$
%
% $head Purpose$$
% This routine solves the following linear system of equations for
% $latex e$$:
% $latex \[
%       A * e = r
% \] $$
% where the symmetric block tridiagonal matrix $latex A$$ is defined by
% $latex \[
% A =
% \left( \begin{array}{ccccc}
% b_1    & c_2^\R{T} & 0          & \cdots & 0      \\
% c_2    & b_2       & c_3^\R{T}  &        & \vdots \\
% \vdots &           & \ddots     &        &        \\
% 0      & \cdots    &            & b_N    & c_N
% \end{array} \right)
% \] $$
%
% The routine $cref ckbs_tridiag_solve_b$$ solves the same problem.
% The difference is that this routine runs the reduction starting
% at the first block rather than the last block of the matrix,
% which is less stable for this application.
%
% $head b$$
% The argument $icode b$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       b_k = b(:,:,k)
% \]$$
% and $icode b$$ has size $latex n \times n \times N$$.
%
% $head c$$
% The argument $icode c$$ is a three dimensional array,
% for $latex k = 2 , \ldots , N$$
% $latex \[
%       c_k = c(:,:,k)
% \]$$
% and $icode c$$ has size $latex n \times n \times N$$.
%
% $head r$$
% The argument $icode r$$ is an $latex (n * N) \times m$$ matrix.
%
% $head e$$
% The result $icode e$$ is an $latex (n * N) \times m$$ matrix.
%
% $head lambda$$
% The result $icode lambda$$ is a scalar equal to the
% logarithm of the determinant of $latex A$$.
%
% $head Reference$$
% $icode The marginal likelihood for parameters in a discrete
% Gauss Markov process$$,
% Bradley M. Bell,
% IEEE Transactions on Signal Processing,
% Vol. 48,
% No. 3,
% March 2000.
%
% $head Lemma$$
% The following result is proven as Lemma 6 of the reference above:
% Suppose that for $latex k = 1 , \ldots , N$$
% $latex \[
% \begin{array}{rcl}
%       b_k & = & u_k + q_{k-1}^{-1} + a_k * q_k^{-1} * a_k^\R{T}  \\
%       c_k & = & q_{k-1}^{-1} * a_k^\R{T}
% \end{array}
% \] $$
% where $latex u_k$$ is symmetric positive semi-definite and
% $latex q_k$$ is symmetric positive definite.
% It follows that the algorithm used by $code ckbs_tridiag_solve$$
% is well conditioned and will not try to invert singular matrices.
%
% $children#
%       example/tridiag_solve_ok.m
% #$$
%
% $head Example$$
% The file $cref tridiag_solve_ok.m$$ contains an example and test of
% $code ckbs_tridiag_solve$$.
% It returns true if $code ckbs_tridiag_solve$$ passes the test
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------
function [e, lambda] = ckbs_tridiag_solve(b, c, r)
    %
    % Determine the size of the problem
    n = size(b, 1);
    m = size(r, 2);
    N = size(b, 3);
    %
    % Initialize
    df     = b;
    sf     = zeros(n, m, N);
    blk_n = 1 : n;
    for k = 1 : N
        sf(:, :, k) = r(blk_n, :);
        blk_n      = blk_n + n;
    end
    %
    % Step 3
    for k = 2 : 1 : N
        sk1        = sf(:, :, k-1);
        dk1        = df(:, :, k-1);
        ck         = c(:, :, k);
        df(:, :, k) = b(:, :, k) - ck * chol_solve(dk1, ck');
        sf(:, :, k) = sf(:, :, k) - ck * chol_solve(dk1, sk1);
    end
    %
    % Step 4
    sf(:, :, N) = chol_solve(df(:,:, N),  sf(:,:, N));
    lambda     = log( det( df(:, :, N) ) );
    %
    % Step 5
    for k = N-1 : -1 : 1
        ek1        = sf(:, :, k+1);
        ck1        = c(:, :, k+1);
        dk         = df(:, :, k);
        sf(:, :, k) = chol_solve(dk,  sf(:, :, k) - ck1' * ek1);
        lambda     = lambda + log( det(dk) );
    end
    %
    % return value
    e     = zeros(n * N, m);
    blk_n = 1 : n;
    for k = 1 : N
        e(blk_n, :) = sf(:, :, k);
        blk_n      = blk_n + n;
    end
    return
end

function [x] = chol_solve(A, y)
    A_chol = chol(A);
    x      = A_chol \ ( A_chol' \ y );
    return
end
