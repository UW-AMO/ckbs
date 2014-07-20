% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin direct_h.m$$ $newlinech %$$
%
% $spell
%       end end
%       Jacobian
%       nl
%       ckbs
%       hk
%       xk
%       sqrt
%       dr
%       params
% $$
%
% $section ckbs_nonlinear: Example Direct Measurement Model$$
% $index measure, direct$$
% $index direct, measure$$
% $index ckbs_nonlinear, measure$$
%
% $head Syntax$$
% $codei%[%hk%] %  %  = direct_h(%k%, %xk%, %params%)
% %$$
% $codei%[%hk%, %Hk%] = direct_h(%k%, %xk%, %params%)
% %$$
%
% $head Notation$$
% $codei%
%       %index%     = params.direct_h_index
%       %m%         = size(%index%, 1)
%       %n%         = size(%xk%, 1)
% %$$
%
% $head Purpose$$
% For $icode%i% = 1 , %...% , %m%$$,
% the mean of the $th i$$ measurement given $icode xk$$
% (the state at time index $icode k$$) is
% $codei%
%       %xk%( %index%(%i%) )
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
% is an integer column vector that specifies which components
% of $icode xk$$ correspond to the direct measurements.
% Each element of $icode index$$ must be between one and $icode n$$.
%
% $head hk$$
% The return value $icode hk$$ is a column vector of length $icode m$$
% with $th i$$ element equal to the mean of the
% $th i$$ measurement given $icode xk$$.
%
% $head Hk$$
% The return value $icode Hk$$ is a $icode%m% x %n%$$ matrix equal to the
% Jacobian of $icode hk$$ w.r.t $icode xk$$.
%
% $head Source Code$$
% $newlinech $$ $codep
function [hk, Hk] = direct_h(k, xk, params)
    index = params.direct_h_index;
    m        = size(index, 1);
    n        = size(xk, 1);
    if (size(xk, 2)~=1) | (size(index,2)~=1)
        size_xk_2    = size(xk, 2)
        size_index_2 = size(index, 2)
        error('direct_h: xk or index is not a column vector')
    end
    if (max(index) > n) | (min(index) < 1)
        max_index = max(index)
        min_index = min(index)
        error('direct_h: max(index) > size(xk, 1) or min(index) < 1')
    end
    hk  = xk(index);
    Hk  = zeros(m, n);
    for i = 1 : m
        Hk( i, index(i) ) = 1;
    end
    return
end
% $$ $newlinech %$$
% $end
