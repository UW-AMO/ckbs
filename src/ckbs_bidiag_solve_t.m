% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_bidiag_solve_t$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_bidiag_solve_t
%       Bidiagonal
% $$
%
% $section  Block Bidiagonal Algorithm$$
%
% $index ckbs_bidiag_solve_t$$
% $index bidiag_solve_t$$
%
% $index solve, block bidiagonal transpose$$
% $index bidiagonal transpose, block solve$$
%
% $head Syntax$$
% $codei/[/e/, /lambda/] = ckbs_bidiag_solve_t(/b/, /c/, /r/)/$$
%
% $head Purpose$$
% This routine solves the following linear system of equations for
% $latex e$$:
% $latex \[
%       A^T * e = r
% \] $$
% where the bidiagonal matrix $latex A$$ is defined by
% $latex \[
% A =
% \left( \begin{array}{ccccc}
% b_1    & 0 & 0          & \cdots       \\
% c_2    & b_2       & 0  &         \\
% \vdots &           & \ddots     &               \\
% 0      & \cdots    & c_N           & b_N    
% \end{array} \right)
% \] $$
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
%       example/bidiag_solve_t_ok.m
% #$$
%
% $head Example$$
% The file $cref bidiag_solve_t_ok.m$$ contains an example and test of
% $code ckbs_bidiag_solve_t$$.
% It returns true if $code ckbs_bidiag_solve_t$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------------
% --
function [e, lambda] = ckbs_bidiag_solve_t(c, a, r)
    %
    % Determine the size of the problem
    n = size(c, 1);
    m = size(r, 2);
    N = size(c, 3);
    %


    % Initialize
    s        = zeros(n, m, N);
    d        = zeros(n, n, N);

    % return value
    e     = zeros(n * N, m);
    blk_n = (N-1)*n + 1:n*N;
    ck = c(:,:,N);
    rk = r(blk_n, :);
    e(blk_n, :) = chol_solve(ck', rk);
    lambda = log(det(ck));

    for k = N-1:-1:1
        ek1 = e(blk_n, :);
        blk_n = blk_n - n;
        ck = c(:,:,k);
        lambda = lambda + log(det(ck));
        rk = r(blk_n, :);
        ak = a(:,:,k+1);
        tempk = rk - ak'*ek1;
        e(blk_n, :) = chol_solve(ck', tempk);
    end

    return
end

function [x] = chol_solve(A, y)
    %eig(A)
    A_chol = chol(A);
    x      = A_chol \ ( A_chol' \ y );
    return
end
