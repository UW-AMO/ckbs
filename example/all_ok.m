% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-11
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin all_ok.m$$  $newlinech %$$
% $spell
%       end end
%       blkbidiag
%       symm
%       pcg
%       nl
%       vanderpol_sim
%       nargin
%       ckbs
%       blkdiag
%       blktridiag
%       Mul
%       ok ok
%       ckbs_sumsq_hes
%       Obj
%       ckbs_tridiag_solve
%       ckbs_tridiag_solve_b
%       ckbs_tridiag_solve_cg
%       invdiag
%       feval
%       vec
%       cputime
%       cputime
%       disp
%       num2str
%       ind
%       findstr
%       Kuhn
%       dotdot
%       bidiag
% $$
%
% $section Run All Correctness Tests$$
%
% $index all_ok$$
% $index test, run all$$
% $index correct, test all$$
%
% $head Syntax$$
% $codei%all_ok%$$
%
% $children%
%       example/test_path.m
% %$$
% $head Test Utility Functions$$
% $table
% $rref test_path.m$$
% $tend
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = all_ok()
    test_path;
    ok = true;
    t0 = cputime;
    ok = ok & one_ok('t_obj_ok');
    ok = ok & one_ok('t_grad_ok');
    ok = ok & one_ok('affine_ok_box');
    ok = ok & one_ok('bidiag_solve_ok');
    ok = ok & one_ok('bidiag_solve_t_ok');
    ok = ok & one_ok('blkbidiag_symm_mul_ok');
    ok = ok & one_ok('blkbidiag_mul_ok');
    ok = ok & one_ok('blkbidiag_mul_t_ok');
    ok = ok & one_ok('t_general_ok');
    ok = ok & one_ok('L1_affine_ok');
    ok = ok & one_ok('L1_nonlinear_ok');
    ok = ok & one_ok('blkdiag_mul_ok');
    ok = ok & one_ok('blkdiag_mul_t_ok');
    ok = ok & one_ok('blktridiag_mul_ok');
    ok = ok & one_ok('kuhn_tucker_ok');
    ok = ok & one_ok('kuhn_tucker_L1_ok');
    ok = ok & one_ok('newton_step_ok');
    ok = ok & one_ok('newton_step_L1_ok');
    ok = ok & one_ok('get_started_ok');
    ok = ok & one_ok('sumsq_grad_ok');
    ok = ok & one_ok('process_grad_ok');
    ok = ok & one_ok('sumsq_hes_ok');
    ok = ok & one_ok('process_hes_ok');
    ok = ok & one_ok('sumsq_obj_ok');
    ok = ok & one_ok('L2L1_obj_ok');
    ok = ok & one_ok('tridiag_solve_ok');
    ok = ok & one_ok('tridiag_inv_ok');
    ok = ok & one_ok('tridiag_solve_b_ok');
    ok = ok & one_ok('tridiag_solve_pcg_ok');
    ok = ok & one_ok('vanderpol_sim_ok');
    ok = ok & one_ok('vanderpol_ok');
    constraint={'no_constraint', 'box_constraint', 'sine_constraint'};
    for c = 1 : 3
        if sine_wave_ok(constraint{c}, false)
            ['Ok:    sine_wave_ok ', constraint{c}]
        else
            ['Error: sine_wave_ok ', constraint{c}]
            ok = false
        end
    end
    %
    t1 = cputime;
    if ok
        disp(['All tests passed: cputime (secs) = ', num2str(t1-t0)]);
    else
        disp(['One or more tests failed: cputime (secs) = ', num2str(t1-t0)]);
    end
    return
    end
    function [ok] = one_ok(name)
    ok = feval(name);
    if ok
        ['Ok:    ', name ]
    else
        ['Error: ', name ]
    end
    return
end
% $$ $newlinech %$$
% $end
