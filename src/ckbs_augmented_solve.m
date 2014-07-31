% Augmented solve: solves \tilde C [x; y] = [r1; r2], where 
% \tilde C = [C B \\ B' D], with 
% C block tridiagonal, having blocks of b on the diagonal, 
% and blocks of c on off-diagonals. 
% B and D are passed in as full matrices. 
function [e1 e2 var] = ckbs_augmented_solve(b, c, B, D, r1, r2)
    
f2 = r2-B'*(ckbs_tridiag_solve(b, c, r1));
    DmBB = D - B'*(ckbs_tridiag_solve(b, c, B));
    
%    fprintf('Covariance estimate is: %7.1e\n', 1/DmBB);
    
    
 %   condNum = cond(DmBB);
  %  fprintf('Condition number of Schur complement is %7.1e\n', condNum);
    
    
    e2 = DmBB\f2;
    e1 = ckbs_tridiag_solve(b, c, r1-B*e2);
    var = DmBB\eye(size(DmBB));
    
%    e = [e1; e2];    
end