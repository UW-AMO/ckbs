% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: sasha.aravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_blkbidiag_mul_t$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       blkbidiag
%       Mul
%       Blk
%       Bdia
%       Boffdia
%       bidiagonal
%       bidiag
% $$
%
% $section Packed Lower Block Bidiagonal Matrix Transpose Times a Vector$$
%
% $index ckbs_blkbidiag_mul$$
% $index blkbidiag_mul$$
%
% $index block, lower bidiagonal transpose multiply$$
% $index lower bidiagonal, block transpose multiply$$
% $index multiply, block lower bidiagonal transpose$$
%
% $head Syntax$$
% $codei%[%w%] = ckbs_blkbidiag_mul_t(%Bdia%, %Boffdia%, %v%)%$$
%
% $head Purpose$$
% This routine enables one to used the packed form of a lower block bidiagonal matrix
% and returns the transpose of the matrix times a vector or matrix; i.e.,
% $latex \[
%       W = B^T * V
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
% $head Boffdia$$
% The argument $icode Boffdia$$ is an $latex m \times n \times N$$ array.
% For $latex k = 2 , \ldots , N$$ we define
% $latex C_k \in \B{R}^{m \times n}$$ by
% $latex \[
%       C_k = Boffdia(:, :, k)
% \] $$
%
%
% $head B$$
% The matrix $latex B$$ is defined by
% $latex \[
% B
% =
% \left( \begin{array}{cccc}
%       B_1 & 0     & 0      &           \\
%       C_2^T   & 0   & \ddots      & 0         \\
%       0   & \ddots      & \ddots & 0         \\
%           & 0      & C_N^T      & B_N
% \end{array} \right)
% \] $$
%
% $head V$$
% The argument $icode V$$ is a matrix of size $latex n * N \times
% k$$, with no restriction on $latex k$$.
%
% $head W$$
% The result $icode W$$ is a matrix of size $latex m * N \times k$$.
%
%
% $children#
%       example/blkbidiag_mul_t_ok.m
% #$$
%
% $head Example$$
% The file $cref blkbidiag_mul_t_ok.m$$ contains an example and test of
% $code ckbs_blkbidiag_mul_t$$.
% It returns true if $code ckbs_blkbidiag_mul_t$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------

% can multiply by n*N x n matrix

function [w] = ckbs_blkbidiag_mul_t(Bdia, Boffdia, v)
    n     = size(Bdia, 1);
    N     = size(Bdia, 3);
    m     = size(v, 2);


    blk_n_new = 1 : n;
    %
    w     = zeros(n * N, m);
    vCur    = zeros(n, m);
    v = [v; zeros(n,m)];

    BNew = zeros(n, n);
    Blast = zeros(n,n);
    Boffdia(:,:,N+1) = zeros(n,n);

    for k = 1 : N;
        blk_n = blk_n_new;
        blk_n_new = blk_n + n;

        vCur = v(blk_n, :);
        vPlus = v(blk_n_new, :);

        BNew = Boffdia(:, :, k+1);

        w(blk_n, :) = Bdia(:,:, k)' * vCur + BNew' * vPlus;
    end
    return
end
