% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-11
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin blkbidiag_symm_mul_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       Bdiag
%       Boffdiag
%       blk
%       blkdiag
%       blktridiag
%       blkbidiag
%       mul
%       bidiag
%       tridiag
%       ak
%       Hdia
%       Hlow
%       Hcheck
%       speye
%       symm
%       Ddia
%       dk
%       Hfull
%       Dfull
%       Ddia
% $$
%
% $section blkbidiag_symm_mul Example and Test$$
%
% $index ckbs_blkbidiag_symm_mul, example and test$$
% $index blkbidiag_symm_mul, example and test$$
% $index example, blkbidiag_symm_mul$$
% $index test, blkbidiag_symm_mul$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = blkbidiag_symm_mul_ok()
ok = true;
% -------------------------------------------------------------
% You can change these parameters
m    = 2;
n    = 3;
N    = 4;
% -------------------------------------------------------------
% Define the problem
rand('seed', 123);

rand('seed', 123);

Hdia  = rand(n, n, N);
Hlow  = rand(n, n, N);
Ddia  = zeros(n, n, N); 

for k = 1:N
    dk = rand(n, n);
    Ddia(:,:,k)  = dk'*dk;
end

Hfull = ckbs_blkbidiag_mul(Hdia, Hlow, eye(n*N));
Dfull = ckbs_blkdiag_mul(Ddia, eye(n*N));

HDHfull = Hfull * Dfull *  Hfull';

[r, s]  = ckbs_blkbidiag_symm_mul(Hdia, Hlow, Ddia);

HDH     = ckbs_blktridiag_mul(r, s, eye(n*N));

check = HDH - HDHfull;

ok    = ok & ( max(max(abs(check))) < 1e-10 );
return
end
% $$ $newlinech %$$
% $end
