% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_blkdiag_mul_t$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       Mul
%       Blk
%       Bdia
% $$
%
% $section Transpose of Packed Block Diagonal Matrix Times a Vector
% or Matrix$$
%
% $index ckbs_blkdiag_mul_t$$
% $index blkdiag_mul_t$$
%
% $index block, diagonal multiply$$
% $index diagonal, block multiply$$
% $index multiply, block diagonal$$
%
% $head Syntax$$
% $codei%[%w%] = ckbs_blkdiag_mul_t(%Bdia%, %v%)%$$
%
% $head Purpose$$
% This routine enables one to used the packed form of a block diagonal matrix
% and returns the transpose of the matrix times a matrix; i.e.,
% $latex \[
%       W = B^\R{T} * V
% \] $$
%
% $head Bdia$$
% The argument $icode Bdia$$ is an $latex m \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex B_k \in \B{R}^{m \times n}$$ by
% $latex \[
%       B_k = Bdia(:, :, k)
% \] $$
%
% $head B$$
% The matrix $latex B$$ is defined by
% $latex \[
% B
% =
% \left( \begin{array}{cccc}
%       B_1 & 0      & 0      &           \\
%       0   & B_2    & 0      & 0         \\
%       0   & 0      & \ddots & 0         \\
%           & 0      & 0      & B_N
% \end{array} \right)
% \] $$
%
% $head V$$
% The argument $icode V$$ is a matrix of size $latex m * N \times
% p$$, with no restrictions on $latex p$$.
%
% $head W$$
% The result $icode W$$ is a matrix of size $latex n * N \times p$$.
%
% $children#
%       example/blkdiag_mul_t_ok.m
% #$$
%
% $head Example$$
% The file $cref blkdiag_mul_t_ok.m$$ contains an example and test of
% $code ckbs_blkdiag_mul_t$$.
% It returns true if $code ckbs_blkdiag_mul_t$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------
function [w] = ckbs_blkdiag_mul_t(Bdia, v)
    m     = size(Bdia, 1);
    n     = size(Bdia, 2);
    N     = size(Bdia, 3);
    p     = size(v, 2);

    blk_m = 1 : m;
    blk_n = 1 : n;
    %
    w     = zeros(n * N, p);
    for k = 1 : N;
        w(blk_n, :) = Bdia(:,:, k)' * v(blk_m, :);
        blk_m = blk_m + m;
        blk_n = blk_n + n;
    end
    return
end
