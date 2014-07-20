% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin no_f.m$$ $newlinech %$$
%
% $spell
%       Jacobian
%       ckbs
%       fk
%       xk
%       nof
% $$
%
% $section ckbs_nonlinear: Example of No Constraint$$
% $index unconstrained, utility$$
% $index constraint, none$$
% $index ckbs_nonlinear, unconstrained$$
%
% $head Syntax$$
% $codei%[%fk%] %  %  = no_f(%k%, %xk%)
% %$$
% $codei%[%fk%, %Fk%] = no_f(%k%, %xk%)
% %$$
%
% $head Purpose$$
% Implements a constraint that is always satisfies by defining the function
% $latex \[
%       f_k ( x_k ) \equiv -1
% \] $$
% so that the condition $latex f_k ( x_k ) \leq 0$$ is satisfied for all
% $latex x_k$$.
%
% $head k$$
% is an integer scalar specifying the time index (not used).
%
% $head xk$$
% The column vector $icode xk$$ specifies the current state vector
% (only used to determine the size of the state vector at each time point).
%
% $head fk$$
% The return value $icode fk$$ is a scalar equal to minus one.
%
% $head Fk$$
% The return value $icode Fk$$ is a row vector equal to the Jacobian of
% $icode fk$$ w.r.t $icode xk$$; i.e. zero.
%
% $head Source Code$$
% $newlinech $$ $codep
function [fk, Fk] = no_f(k, xk)
    % no constraint case: f(x) = -1 for all x
    n  = size(xk, 1);
    fk = -1;
    Fk = zeros(1, n);
    return
end
% $$ $newlinech %$$
% $end
