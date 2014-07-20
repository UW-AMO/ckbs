% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: sasha.aravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_blkbidiag_symm_mul$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       blkbidiag
%       Mul
%       Blk
%       Bdia
%       Boffdia
%       bidiag
%       blkbidiag
%       tridiagonal
%       tridiag
%       symm
%       Ddia
%       dk
%       Hfull
%       Dfull
%       Ddia
%       Bidiagonal
% $$
%
% $section Packed Block Bidiagonal Matrix Symmetric Multiply$$
%
% $index ckbs_blkbidiag_symm_mul$$
% $index blkbidiag_symm_mul$$
%
% $index block, bidiagonal symmetric multiply$$
% $index bidiagonal symmetric, block multiply$$
% $index symmetric multiply, block bidiagonal$$
%
% $head Syntax$$
% $codei%[%r%, %s%] = ckbs_blkbidiag_symm_mul(%Bdia%, %Boffdia%, %Ddia%)%$$
%
% $head Purpose$$
% This routine takes a packed block bidiagonal matrix and
% a packed block diagonal matrix (optional), and returns
% the symmetric product matrix (as packed tridiagonal matrix):
% $latex \[
%       W = B * D * B^T
% \] $$
% The actual output consists of a packed representation of the
% diagonal and off-diagonal blocks of the matrix $latex W$$.
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
% $head Ddia$$
% The argument $icode Ddia$$ is an $latex n \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex D_k \in \B{R}^{n \times n}$$ by
% $latex \[
%       D_k = Ddia(:, :, k)
% \] $$
%
%
% $head B$$
% The matrix $latex B$$ is defined by
% $latex \[
% B
% =
% \left( \begin{array}{cccc}
%       B_1 & 0      & 0      &           \\
%       C_2^T   & B_2    & \ddots      & 0         \\
%       0   & \ddots      & \ddots & 0         \\
%           & 0      & C_N^T      & B_N
% \end{array} \right)
% \] $$
%
% $head D$$
% The matrix $latex D$$ is defined by
% $latex \[
% D
% =
% \left( \begin{array}{cccc}
%       D_1 & 0      & 0      &           \\
%       0   & D_2    & \ddots      & 0         \\
%       0   & \ddots      & \ddots & 0         \\
%           & 0      & 0      & D_N
% \end{array} \right)
% \] $$
%
%
%
% $head r$$
% The result $icode r$$ is an $latex n \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex r_k \in \B{R}^{n \times n}$$ by
% $latex \[
%       r_k = r(:, :, k)
% \] $$
%
%% $head s$$
% The result $icode s$$ is an $latex m \times n \times N$$ array.
% For $latex k = 2 , \ldots , N$$ we define
% $latex r_k \in \B{R}^{m \times n}$$ by
% $latex \[
%       s_k = s(:, :, k)
% \] $$
%
% $children#
%       example/blkbidiag_symm_mul_ok.m
% #$$
%
% $head Example$$
% The file $cref blkbidiag_symm_mul_ok.m$$ contains an example and test of
% $code ckbs_blkbidiag_symm_mul$$.
% It returns true if $code ckbs_blkbidiag_symm_mul$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------

function [r s] = ckbs_blkbidiag_symm_mul(Bdia, Boffdia, Ddia)

n     = size(Bdia, 1);
m     = size(Boffdia, 1);
N     = size(Bdia, 3);

if nargin < 3
    Ddia = zeros(n, n, N);
    for k = 1:N
        Ddia(:,:,k) = eye(n);
    end
    
end

r     = zeros(n, n, N);
s     = zeros(m, n, N);

bdk    = Bdia(:,:,1);
dk     = Ddia(:,:,1);

r(:,:,1) = Bdia(:,:,1) * Ddia(:,:,1) * Bdia(:,:,1)';

for k = 2 : N
    r(:,:,k) = Boffdia(:,:,k) * Ddia(:,:,k-1) * Boffdia(:,:,k)' ...
        + Bdia(:,:,k) * Ddia(:,:,k) * Bdia(:,:,k)';
    
    s(:,:,k) = Boffdia(:,:,k) * Ddia(:,:,k-1) * Bdia(:,:,k-1)';
end
return
end
