% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin persist_g.m$$ $newlinech %$$
%
% $spell
%       Jacobian
%       ckbs
%       gk
%       xk
%       dt
%       params
% $$
%
% $section ckbs_nonlinear: Example of Persistence Transition Function$$
%
% $index dynamics, persistence$$
% $index transition, persistence$$
% $index persistence, transition$$
% $index ckbs_nonlinear, persistence$$
% $index random walk, transition$$
%
% $head Syntax$$
% $codei%[%gk%] %  %  = persist_g(%k%, %xk1%, %params%)
% %$$
% $codei%[%gk%, %Gk%] = persist_g(%k%, %xk1%, %params%)
% %$$
%
% $head Notation$$
% $codei%
%    %initial% = persist_g.initial
%       %n%       = size(%xk1%, 1)
% %$$
%
% $head Purpose$$
% Implements the persistence model for state transitions; i.e.,
% the mean of the state at time index $icode k$$ given the state at time
% index $icode%k%-1%$$ is its value at time index $icode%k%-1%$$.
% (This corresponds to a random walk model.)
%
% $head initial$$
% is a column vector of length $icode n$$ specifying the initial estimate
% for the state vector at time index one.
%
% $head k$$
% is a positive integer scalar specifies the current time index.
%
% $head xk1$$
% is a column vector specifying a value for
% the state vector at the previous time index $icode%k%-1%$$.
%
% $head gk$$
% If $icode%k% == 1%$$,
% the return value $icode gk$$ is equal to $icode initial$$.
% Otherwise, $icode gk$$ is equal to $icode xk1$$.
%
% $head Gk$$
% The return value $icode Gk$$ is an $icode%n% x %n%$$ matrix equal to the
% Jacobian of $icode gk$$ w.r.t $icode xk1$$; i.e., the identity matrix.
%
% $head Source Code$$
% $newlinech $$ $codep
function [gk, Gk] = persist_g(k, xk1, params)
    initial = params.persist_g_initial;
    n       = size(xk1, 1);
    if k == 1
        gk = initial;
        Gk = zeros(n, n);
    else
        gk = xk1;
        Gk = eye(n);
    end
    return
end
% $$ $newlinech %$$
% $end
