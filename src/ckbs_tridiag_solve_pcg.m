% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2012
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:   saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_tridiag_solve_pcg$$ $newlinech %$$
% $spell
%       ckbs
%       iter
%       Matlab
%       cg
%       blkdiag
%       ckbs_tridiag_solve_pcg
%       Tridiagonal
%       Bradley
% $$
%
% $section Symmetric Block Tridiagonal Algorithm (Conjugate Gradient version)$$
%
% $index ckbs_tridiag_solve_pcg$$
% $index tridiag_solve_pcg$$
%
% $index solve, cg block tridiagonal$$
% $index tridiagonal, pcg block solve$$
% $index cg, block tridiagonal solve$$
%
% $head Syntax$$
% $codei/[/e/, /iter/] = ckbs_tridiag_solve_pcg(/b/, /c/, /r/)/$$
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
% The routines $cref ckbs_tridiag_solve$$ and $cref ckbs_tridiag_solve_b$$
% solve the same problem, but only for one RHS. It is basically a wrapper
% for Matlab's pcg routine. 
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
% $head in$$
% The argument $icode in$$ is an $latex (n * N) \times 1$$ vector.
%
% $head e$$
% The result $icode e$$ is an $latex (n * N) \times 1$$ vector.
%
% $head iter$$
% The result $icode iter$$ is a scalar equal to the
% number of iterations that the cg method took to finish.
%
%
% $children#
%       example/tridiag_solve_pcg_ok.m
% #$$
%
% $head Example$$
% The file $cref tridiag_solve_pcg_ok.m$$ contains an example and test of
% $code ckbs_tridiag_solve_cg$$.
% It returns true if $code ckbs_tridiag_solve_pcg$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------------
% --
function [e, flag, relres, iter] = ckbs_tridiag_solve_pcg(c, a, in)
    %    

    n = size(c,1);
    N = size(c,3);
    [e flag relres iter]     =  pcg(@(x)ckbs_blktridiag_mul(c, a, x), in, [], n*N);
   

    return
end

