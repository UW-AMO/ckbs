% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin sine_f.m$$ $newlinech %$$
%
% $spell
%       cos
%       Jacobian
%       nl
%       ckbs
%       fk
%       xk
%       params
% $$
%
% $section ckbs_nonlinear: Example of Nonlinear Constraint$$
% $index sine wave, constraint$$
% $index constraint, sine wave$$
% $index ckbs_nonlinear, sine wave constraint$$
%
% $head Syntax$$
% $codei%[%fk%] %  %  = sine_f(%k%, %xk%, %params%)
% %$$
% $codei%[%fk%, %Fk%] = sine_f(%k%, %xk%, %params%)
% %$$
%
% $head Notation$$
% $codei%
%       %index%  = params.sine_f_index
%       %offset% = params.sine_f_offset
% %$$
%
% $head Purpose$$
% Implements an upper limit that is an offset sine wave
% To be specific,
% $codei%
% %xk%( %index%(%2%) ) <= sin( %xk%( %index%(1) + %offset%(1) ) ) + %offset%(2)
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
% is an integer column vector with two elements
% specifying the state vector indices for the sine wave constraint
% (see purpose above).
% Each such index must be between one and $icode n$$.
%
% $head offset$$
% is an integer column vector with two elements
% specifying the offsets for the sine wave (see purpose above).
%
% $head fk$$
% The return value $icode fk$$ is a scalar equal to
% $codei%
% %xk%( %index%(%2%) ) - sin( %xk%( %index%(1) + %offset%(1) ) ) - %offset%(2)
% %$$
% so that the condition $icode%fk% <= 0%$$ is equivalent to
% the constraints in the purpose above.
%
% $head Fk$$
% The return value $icode Fk$$ is a row vector with $icode%n%$$ elements
% equal to the derivative of $icode fk$$ w.r.t $icode xk$$.
%
% $head Source Code$$
% $newlinech $$ $codep
function [fk, Fk] = sine_f(k, xk, params)
    n      = size(xk, 1);
    index  = params.sine_f_index;
    offset = params.sine_f_offset;
    if (size(index, 1) ~= 2) | (size(index, 2) ~= 1)
        size_index_1 = size(index, 1)
        size_index_2 = size(index, 2)
        error('sine_f: index is not a column vector with two elements')
    end
    if (size(offset, 1) ~= 2) | (size(offset, 2) ~= 1)
        size_offset_1 = size(offset, 1)
        size_offset_2 = size(offset, 2)
        error('sine_f: offset is not a column vector with two elements')
    end
    if (max(index) > n) | (min(index) < 1)
        max_index = max(index)
        min_index = min(index)
        error('sine_f: max(index) > size(xk, 1) or min(index) < 1')
    end
    fk = xk( index(2) ) - sin( xk( index(1) + offset(1) ) ) - offset(2);
    Fk           = zeros(1, n);
    Fk(index(2)) = 1;
    Fk(index(1)) = Fk(index(1)) - cos( xk( index(1) + offset(1) ) );
    return
end
% $$ $newlinech %$$
% $end
