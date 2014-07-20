% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2014
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    sasha.aravkin at gmail dot com
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_tridiag_inv$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_tridiag_inv
%       Tridiagonal
%       var
%       ifull
%       ind
% $$
%
% $section Symmetric Inverse Covariance Block Tridiagonal Algorithm$$
%
% $index ckbs_tridiag_inv$$
% $index tridiag_inv$$
%
% $index inv, block tridiagonal$$
% $index tridiagonal, block inv$$
%
% $head Syntax$$
% $codei/[/var/] = ckbs_tridiag_inv(/b/, /c/)/$$
%
% $head Purpose$$
% This routine computes the block diagonal of the inverse 
% of a block tridiagonal system $latex A$$ is defined by
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
%
% $head var$$
% The result $icode var$$ is an $latex (n * n * N)$$ array.
%
%
%
% $children#
%       example/tridiag_inv_ok.m
% #$$
%
% $head Example$$
% The file $cref tridiag_inv_ok.m$$ contains an example and test of
% $code ckbs_tridiag_inv$$.
% It returns true if $code ckbs_tridiag_inv$$ passes the test
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------
function [var] = ckbs_tridiag_inv(b, c)
%
% Determine the size of the problem
n = size(b, 1);
N = size(b, 3);
%

% Initialize 
d     = b;

var   = zeros(n, n, N);

for k = 2 : 1 : N
	dk1        = d(:, :, k-1);
	ck         = c(:, :, k);
	d(:, :, k) = b(:, :, k) - ck * chol_solve(dk1, ck');    
end


var(:, :, N) = chol_solve(d(:,:, N),  eye(n));

for k = N-1 : -1 : 1
    dk_inv = chol_solve(d(:,:,k), eye(n));
    dat    = dk_inv * c(:,:,k+1)';
    var(:,:,k) = dk_inv + dat*var(:,:,k+1)*dat';
end
return
end

function [x] = chol_solve(A, y)
%eig(A)
A_chol = chol(A);
x      = A_chol \ ( A_chol' \ y );
return
end
