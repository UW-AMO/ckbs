load 'TridiagBreakEx'

y = ckbs_tridiag_solve_b(D, A, -dVec);

%works!

dVecEst = -ckbs_blktridiag_mul(D, A, y);

norm(dVec - dVecEst);


y = ckbs_tridiag_solve(D, A, -dVec);


% Fails!
