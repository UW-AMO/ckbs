% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin distance_h.m$$ $newlinech %$$ 
%
% $spell
%	end end
%	Jacobian
%	nl
%	ckbs
%	hk
%	xk
%	sqrt
%	dr
%   params
% $$
%
% $section ckbs_nonlinear: Example of Distance Measurement Model$$
% $index measure, distance$$
% $index distance, measure$$
% $index ckbs_nonlinear, measure$$
%
% $head Syntax$$
% $codei%[%hk%] %  %  = distance_h(%k%, %xk%, %params%)
% %$$
% $codei%[%hk%, %Hk%] = distance_h(%k%, %xk%)
% %$$
%
% $head Notation$$
% $codei%
%	%index%     = distance_h.index
%	%position%  = distance_h.position
%	%n%         = size(%xk%, 1)
%	%m%         = size(%position%, 2)
% %$$
%
% $head Purpose$$
% For $icode%i% = 1 , %...% , %m%$$,
% the mean of the $th i$$ measurement given $icode xk$$
% (the state at time index $icode k$$) is
% the distance from $icode xk$$ to $icode%position%(:,%i%)%$$; i.e.
% $codei%
%	norm( %xk%(%index%) - %position%(:, %i%) )
% %$$ 
%
% $head index$$
% the integer column vector $icode index$$ specifies which components
% of $icode xk$$ correspond to the current moving position which the distance
% measurements are made to.
% It must have the same number of rows as the matrix $icode position$$.
%
% $head position$$
% The matrix $icode position$$ with $icode m$$ columns.
% Each column specifies a fixed locations that a
% distance measurements are made from.
%
% $head k$$
% is an integer scalar specifying the time index (not used).
%
% $head xk$$
% is a column vector with length $icode n$$ specifying a value for
% the state vector at the current time index.
%
% $head hk$$
% The return value $icode hk$$ is a column vector of length $icode m$$
% with $th i$$ element equal to the mean of the 
% distance from $icode%position%(:,%i%)%$$ to $icode%xk%(%index%)%$$. 
%
% $head Hk$$
% The return value $icode Hk$$ is a $icode%m% x %n%$$ matrix equal to the
% Jacobian of $icode hk$$ w.r.t $icode xk$$.
% (In the special case where one of the distances is zero, 
% the corresponding row of $icode HK$$ is returned as zero.)
%
% $head Source Code$$
% $newlinech $$ $codep
function [hk, Hk] = distance_h(k, xk, params)
	index    = params.distance_h_index;
	position = params.distance_h_position;
	n        = size(xk, 1);
	m        = size(position, 2);
	if (size(xk, 2)~=1) | (size(index,2)~=1)
		size_xk_2    = size(xk, 2)
		size_index_2 = size(index, 2)
		error('distance_h: xk or index is not a column vector')
	end
	if size(position,1) ~= size(index,1)
		size_position_1 = size(position, 1)
		size_index_1    = size(index, 1)
		error('distance_h: position and index have different number of rows')
	end
	hk  = zeros(m, 1);
	Hk  = zeros(m, n);
	p_x = xk(index);
	for i = 1 : m
		p_i          = position(:, i);
		hk(i)        = norm( p_x - p_i );
		if hk(i) ~= 0
			Hk(i, index) = (p_x - p_i) / hk(i);
		end
	end
	return
end
% $$ $newlinech %$$
% $end
