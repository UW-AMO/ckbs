% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin pos_vel_g.m$$ $newlinech %$$
%
% $spell
%       Jacobian
%       ckbs
%       gk
%       xk
%       dt
%       pos_vel
%       params
% $$
%
% $section ckbs_nonlinear: Example Position and Velocity Transition Model$$
% $index dynamics, position and velocity$$
% $index transition, position and velocity$$
% $index position, velocity transition$$
% $index velocity, position transition$$
% $index ckbs_nonlinear, transition$$
%
% $head Syntax$$
% $codei%[%gk%] %  %  = pos_vel_g(%k%, %xk1%, %params%)
% %$$
% $codei%[%gk%, %Gk%] = pos_vel_g(%k%, %xk1%)
% %$$
%
% $head Notation$$
% $codei%
%    %initial% = params.pos_vel_initial
%       %dt%      = params.pos_vel_g_dt
%       %n%       = size(%xk1%, 1)
%       %v%(%j%)  = %xk1%(2 * %j% - 1)
%       %p%(%j%)  = %xk1%(2 * %j%)
% %$$
%
% $head Purpose$$
% For $icode%j% = 1 , %...% , %n% / 2%$$,
% $icode%p%(%j%)%$$
% is a position estimate at time index $icode%k%-1%$$, and
% $icode%v%(%j%)%$$
% is the corresponding velocity estimate at time index $icode%k%-1%$$,
% The mean for the corresponding position at time index $icode k$$,
% given $icode xk1$$, is
% $codei%
%       %p%(%j%) + %v%(%j%) * %dt%
% %$$.
% The mean for the corresponding velocity at time index $icode k$$,
% given $icode xk1$$, is
% $codei%
%       %v%(%j%)
% %$$.
%
% $head dt$$
% is a scalar specifying the time between points for this Kalman smoother.
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
% Otherwise, $icode gk$$ is a column vector of length $icode n$$
% equal to the mean for the state at time index $icode k$$ given the
% value for the state at time index $icode%k%-1%$$; i.e., $icode xk1$$.
%
% $head Gk$$
% The return value $icode Gk$$ is an $icode%n% x %n%$$ matrix equal to the
% Jacobian of $icode gk$$ w.r.t $icode xk1$$.
%
% $head Source Code$$
% $newlinech $$ $codep
function [gk, Gk] = pos_vel_g(k, xk1, params)
    dt      = params.pos_vel_g_dt;
    initial = params.pos_vel_g_initial;
    n  = size(xk1, 1);
    %
    if (size(xk1, 2)~=1) | (size(initial, 1)~=n) | (size(initial, 2)~=1)
        size_xk1_1     = size(xk1, 1)
        size_initial_1 = size(initial, 1)
        size_initial_2 = size(initial, 2)
        error('pos_vel_g: initial or xk1 not column vectors with same size')
    end
    %
    Gk = zeros(n, n);
    if k == 1
        gk = initial;
        return;
    end
    % gk(2*j-1) = xk1(2*j-1)
    Gk     = eye(n);
    for j2 = 2 : 2 : n
        % gk(2*j) = xk1(2*j) + xk1(2*j-1) * dt
        Gk(j2, j2-1) = 5*dt;
    end
    gk = Gk * xk1;
    return
end
% $$ $newlinech %$$
% $end
