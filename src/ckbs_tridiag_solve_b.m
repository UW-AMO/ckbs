% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2014
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    sasha.aravkin at gmail dot com
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_tridiag_solve_b$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_tridiag_solve_b
%       Tridiagonal
%       Bradley
% $$
%
% $section Symmetric Block Tridiagonal Algorithm (Backward version)$$
%
% $index ckbs_tridiag_solve_b$$
% $index tridiag_solve$$
%
% $index solve, backward block tridiagonal$$
% $index tridiagonal, backward block solve$$
% $index backward, block tridiagonal solve$$
%
% $head Syntax$$
% $codei/[/e/, /lambda/] = ckbs_tridiag_solve_b(/b/, /c/, /r/)/$$
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
% The routine $cref ckbs_tridiag_solve$$ solves the same problem.
% The difference is that this routine runs the reduction starting
% at the last block rather than the first block of the matrix,
% which is more stable for this application.
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
%
% $children#
%       example/tridiag_solve_b_ok.m
% #$$
%
% $head Example$$
% The file $cref tridiag_solve_b_ok.m$$ contains an example and test of
% $code ckbs_tridiag_solve_b$$.
% It returns true if $code ckbs_tridiag_solve_b$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------------
% --
function [e, lambda] = ckbs_tridiag_solve_b(b, c, r)
    %
    % Determine the size of the problem
    n = size(b, 1);
    m = size(r, 2);
    N = size(b, 3);
    %


    % Initialize
    sb        = zeros(n, m, N);
    db        = zeros(n, n, N);


    blk_n = (n*(N-1) + 1):n*N;
    rk = r(blk_n, :);
    sk = rk;
    sb(:,:,N) = sk;
    dk = b(:,:,N);
    db(:,:,N) = dk;
    ck = c(:,:,N);

    for k = N-1 : -1 : 1
        ck1 = ck;
        ck = c(:,:,k);
        dk1 = dk;
        sk1 = sk;
        blk_n = blk_n - n;
        rk = r(blk_n, :);
        bk = b(:,:,k);
        %  'first chol_solve'
        dk = bk - ck1'* chol_solve(dk1, ck1);
        db(:,:,k) = dk;
        %   'second chol_solve'
        sk = rk - ck1'* chol_solve(dk1, sk1);
        sb(:,:,k) = sk;
    end

    % return value
    e     = zeros(n * N, m);
    blk_n = 1:n;
    dk = db(:,:,1);
    sk = sb(:,:,1);
    %'third chol_solve'
    e(blk_n, :) = chol_solve(dk, sk);
    lambda = log(det(dk));

    for k = 2:1:N
        ek1 = e(blk_n, :);
        blk_n = blk_n + n;
        dk = db(:,:,k);
        lambda = lambda + log(det(dk));
        sk = sb(:,:,k);
        ck = c(:,:,k);
        tempk = sk - ck*ek1;
        %    'fourth chol_solve'
        e(blk_n, :) = chol_solve(dk, tempk);
    end

    return
end

function [x] = chol_solve(A, y)
    %eig(A)
    A_chol = chol(A);
    x      = A_chol \ ( A_chol' \ y );
    return
end
