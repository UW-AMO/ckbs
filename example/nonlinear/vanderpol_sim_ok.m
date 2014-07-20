% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin vanderpol_sim_ok.m$$ $newlinech %$$
%
% $spell
%       vanderpol_sim
%       mu
%       xi
% $$
%
% $section Example Use of vanderpol_sim$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = vanderpol_sim_ok()
    mu     = 1.;
    xi     = [ 2 ; 0 ];
    n_out  = 21;
    step   = .1;
    x_out  = vanderpol_sim(mu, xi, n_out, step);
    t_sol  = [ 0.0    , 0.5    , 1.0    ,  1.5   , 2.0    ];
    x1_sol = [ 2.0000 , 1.8377 , 1.5081 ,  1.0409, 0.3233 ];
    ok     = 1;
    for k = 1 : 5
        ok = ok & abs( x_out(1, (k-1)*5 + 1) - x1_sol(k) ) < 1e-4;
    end
    return
end
% $$ $newlinech %$$
% $end
