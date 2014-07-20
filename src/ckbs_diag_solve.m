% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_diag_solve$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_diag_solve
%       diagonal
% $$
%
% $section  Block Diagonal Algorithm$$
%
% $index ckbs_diag_solve$$
% $index diag_solve$$
%
% $index solve, block diagonal$$
% $index diagonal, block solve$$
%
% $head Syntax$$
% $codei/[/e/, /lambda/] = ckbs_diag_solve(/b/, /r/)/$$
%
% $head Purpose$$
% This routine solves the following linear system of equations for
% $latex e$$:
% $latex \[
%       A * e = r
% \] $$
% where the diagonal matrix $latex A$$ is defined by
% $latex \[
% A =
% \left( \begin{array}{ccccc}
% b_1    & 0 & 0          & \cdots       \\
% 0    & b_2       & 0  &         \\
% \vdots &           & \ddots     &               \\
% 0      & \cdots    & 0           & b_N    
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
%       example/diag_solve_ok.m
% #$$
%
% $head Example$$
% The file $cref diag_solve_ok.m$$ contains an example and test of
% $code ckbs_diag_solve$$.
% It returns true if $code ckbs_diag_solve$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------------
% --
function [e, lambda] = ckbs_diag_solve(c, r)
    %
    % Determine the size of the problem
    n = size(c, 1);
    m = size(r, 2);
    N = size(c, 3);
    %


    % Initialize

    % return value
    e     = zeros(n * N, m);
    blk_n = 1:n;
    ck = c(:,:,1);
    rk = r(blk_n, :);
    e(blk_n, :) = chol_solve(ck, rk);
    lambda = log(det(ck));

    for k = 2:1:N
        blk_n = blk_n + n;
        ck = c(:,:,k);
        lambda = lambda + log(det(ck));
         rk = r(blk_n, :);
        e(blk_n, :) = chol_solve(ck, rk);
    end

    return
end

function [x] = chol_solve(A, y)
    %eig(A)
    A_chol = chol(A);
    x      = A_chol \ ( A_chol' \ y );
    return
end
