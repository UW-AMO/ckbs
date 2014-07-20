% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin box_f.m$$ $newlinech %$$
%
% $spell
%       Jacobian
%       nl
%       ckbs
%       fk
%       xk
%       params
% $$
%
% $section ckbs_nonlinear: Example of Box Constraints$$
% $index box, constraint utility$$
% $index constraint, box utility$$
% $index ckbs_nonlinear, box constraint$$
%
% $head Syntax$$
% $codei%[%fk%] %  %  = box_f(%k%, %xk%, %params%)
% %$$
% $codei%[%fk%, %Fk%] = box_f(%k%, %xk%, %params%)
% %$$
%
% $head Notation$$
% $codei%
%       %lower% = params.box_f_lower
%       %upper% = params.box_f_upper
%       %index% = params.box_f_index
%       %ell%   = 2 * size(%index%, 1)
% %$$
%
% $head Purpose$$
% Implements box constraints with upper and lower limits that are
% the same for all time indices $icode k$$.
% To be specific, for $icode%p% = 1 , %...% , %ell% / 2%$$,
% $codei%
%       %lower%(%p%) <= %xk%( %index%(%p%) ) <= %upper%(%p%)
% %$$
%
% $head k$$
% is an integer scalar specifying the time index (not used).
%
% $head xk$$
% is a column vector with length $icode n$$ specifying a value for
% the state vector at the current time index.
%
% $head index$$
% is an integer column vector with $icode%ell% / 2%$$ elements
% specifying the state vector indices for which
% there is a box constraints.
% Each such index must be between one and $icode n$$.
%
% $head lower$$
% is a column vector with $icode%ell% / 2%$$ elements
% specifying the lower limits for the box constraints.
%
% $head upper$$
% is a column vector with $icode%ell% / 2%$$ elements
% specifying the upper limits for the box constraints.
%
% $head fk$$
% The return value $icode fk$$ is a column vector of length $icode ell$$
% such that the condition $codei%all( %fk% <= 0%$$ is equivalent to
% the constraints in the $cref/purpose/box_f.m/Purpose/$$ above.
%
% $head Fk$$
% The return value $icode Fk$$ is an $icode%ell% x %n%$$ matrix equal to the
% Jacobian of $icode fk$$ w.r.t $icode xk$$.
%
% $head Source Code$$
% $newlinech $$ $codep
function [fk, Fk] = box_f(k, xk, params)
    index = params.box_f_index;
    lower = params.box_f_lower;
    upper = params.box_f_upper;
    ell    = 2 * size(lower, 1);
    n      = size(xk,    1);
    %
    if (size(lower, 2) ~= 1) | (size(upper, 2) ~= 1) | ...
            (size(xk, 2) ~= 1)
        size_lower_2 = size(lower, 2)
        size_upper_2 = size(upper, 2)
        size_index_2 = size(index, 2)
        size_xk_2    = size(xk,    2)
        error('box_f: either lower, upper, index or xk not a column vector')
    end
    if (2 * size(upper, 1) ~= ell) | (2 * size(index, 1) ~= ell)
        size_lower_1 = size(lower, 1)
        size_upper_1 = size(upper, 1)
        size_index_1 = size(index, 1)
        error('box_f: lower, upper, or index has a different size')
    end
    if (max(index) > n) | (min(index) < 1)
        max_index = max(index)
        min_index = min(index)
        error('box_f: max(index) > size(xk, 1) or min(index) < 1')
    end
    %
    fk = zeros(ell, 1);
    Fk = zeros(ell, n);
    for p = 1 : (ell / 2)
        j                = index(p);
        fk(2 * p - 1, 1) = xk(j) - upper(p);
        Fk(2 * p - 1, j) = +1.;
        %
        fk(2 * p,     1) = lower(p) - xk(j);
        Fk(2 * p,     j) = -1.;
    end
    return
end
% $$ $newlinech %$$
% $end
