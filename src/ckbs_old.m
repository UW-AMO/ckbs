% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs$$ $newlinech %$$
% $spell
%	obj
%	ckbs
%	qinv
%	rinv
%	feval
%	gk
%	hk
%	xk
% $$
% $latex \newcommand{\T}{{\rm T}}$$
% $latex \newcommand{\R}{{\bf R}}$$
%
% $section The Constrained Kalman-Bucy Smoother$$
%
% $head Syntax$$
% $codei/[/x_out/, /u_out/, /info/] = ckbs(/
%	g_fun/, /h_fun/, /epsilon/, /x_in/, /z/, /b/, /db/, /qinv/, /rinv/)/$$
%
% $head Purpose$$
% This routine minimizes the
% Kalman-Bucy smoother residual sum of squares objective 
% subject to an affine inequality constraint.
%
% $head Notation$$
% The Kalman-Bucy smoother residual sum of squares is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = & 
% \frac{1}{2} 
% [ z_k - h_k ( x_k ) ]^\T * R_k^{-1} * [ z_k - h_k ( x_k ) ]
% \\
% & + &
% \frac{1}{2} 
% [ x_k - g_k ( x_{k-1} ) ]^\T * Q_k^{-1} * [ x_k - g_k ( x_{k-1} ) ]
% \\
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are 
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
% Note that $latex g_1 ( x_0 )$$ is the initial state estimate
% and $latex Q_1$$ is the corresponding covariance.
%
% $head Problem$$
% Our constrained Kalman-Bucy smoother problem is
% $latex \[
% \begin{array}{rll}
% {\rm minimize} 
%	& S( x_1 , \ldots , x_N ) 
%	& {\rm w.r.t.} \; x_1 \in \R^n , \ldots , x_N \in \R^n
% \\
% {\rm subject \; to} 
%	& b_k + B_k * x_k  \leq 0 
%	& {\rm for} \; k = 1 , \ldots , N
% \end{array}
% \] $$
%
% $head First Order Conditions$$
% A state sequence $latex ( x_1 , \ldots , x_N )$$ is considered a solution
% if there is a Lagrange multiplier sequence $latex ( u_1 , \ldots , u_N )$$ 
% such that the following conditions are satisfied.
% $latex \[
% \begin{array}{rcl}
%	b_k + B_k * x_k                   & \leq & 0           \\
%	0                                 & \leq & u_k         \\
%       | ( B_k^T * u_k + d_k )_j |       & \leq & \varepsilon \\  
%       (u_k)_i * (- b_k - B_k * x_k)_i   & \leq & \varepsilon
% \end{array}
% \] $$
% for $latex j = 1 , \ldots , n$$ and $latex i = 1 , \ldots , \ell$$.
% Here 
% $latex d_k$$ is the partial derivative of $latex S ( x_1 , \ldots , x_N )$$
% with respect to $latex x_k$$
% and $latex (u_k)_i$$ denotes the $th i$$ component of $latex u_k$$.
%
% $head epsilon$$
% The positive scalar $icode epsilon$$ specifies the convergence 
% criteria value; i.e.,
% $latex \[
%	\varepsilon = epsilon
% \] $$
% 
% $head h_fun$$
% The argument $icode h_fun$$ is a function that supports both of the
% following syntaxes
% $codei/
%	[/hk/] = feval(/h_fun/, /k/, /xk/)
%	[/hk/, /Hk/] = feval(/h_fun/, /k/, /xk/)
% /$$
%
% $subhead k$$
% The $icode h_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$.
%
% $subhead xk$$
% The $icode h_fun$$ argument $icode xk$$ is an column vector with
% length $latex n$$. It specifies a value for $latex x_k$$.
%
% $subhead hk$$
% The $icode h_fun$$ result $icode hk$$ is an column vector with
% length $latex m$$ and
% $latex \[
%	hk = h_k ( x\_k ) 
% \] $$
%
% $subhead Hk$$
% If the $icode h_fun$$ result $icode Hk$$ is present in the syntax,
% it is the $latex m \times n$$ matrix given by
% and
% $latex \[
%	Hk = h_k^{(1)} ( x\_k ) 
% \] $$
% 
% $head g_fun$$
% The argument $icode g_fun$$ is a function that supports both of 
% the following syntaxes
% $codei/
%	[/gk/] = feval(/g_fun/, /k/, /xk1/)
%	[/gk/, /Gk/] = feval(/g_fun/, /k/, /xk1/)
% /$$
%
% $subhead k$$
% The $icode g_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$. 
%
% $subhead xk1$$
% The $icode g_fun$$ argument $icode xk1$$ is an column vector with
% length $latex n$$. 
% It specifies a value for $latex x_{k-1}$$.
%
% $subhead gk$$
% The $icode g_fun$$ result $icode gk$$ is an column vector with
% length $latex n$$ and
% $latex \[
%	gk = g_k ( x_{k-1} ) 
% \] $$
%
% $subhead Gk$$
% If the $icode g_fun$$ result $icode Gk$$ is present in the syntax,
% it is the $latex n \times n$$ matrix given by
% and
% $latex \[
%	Gk = g_k^{(1)} ( x_{k-1} ) 
% \] $$
% 
% $head x_in$$
% The argument $icode x_in$$ contains a strictly feasible sequence; 
% to be specific
% for $latex k = 1 , \ldots , N$$
% $latex \[
% 	b(:, k) + db(:,:, k) * x\_in(:, k) < 0
% \]$$
% and $icode x_in$$ is a two dimensional array with size $latex n \times N$$.
% 
% $head z$$
% The argument $icode z$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
% 	z_k = z(:, k)
% \]$$
% and $icode z$$ has size $latex m \times N$$.
% 
% $head b$$
% The argument $icode b$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
% 	b_k = b(:, k)
% \]$$
% and $icode b$$ has size $latex \ell \times N$$.
% If $latex \ell = 0$$, the problem is not constrained; i.e.,
% it is the affine Kalman-Bucy smoother problem.
% 
% $head db$$
% The argument $icode db$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
% 	B_k = db(:,:,k)
% \]$$
% and $icode db$$ has size $latex \ell \times n \times N$$.
% 
% $head qinv$$
% The argument $icode qinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
% 	Q_k^{-1} = qinv(:,:, k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
% 
% $head rinv$$
% The argument $icode rinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
% 	R_k^{-1} = rinv(:,:, k)
% \]$$
% and $icode rinv$$ has size $latex m \times m \times N$$.
% 
% $head x_out$$
% The result $icode x_out$$ contains the optimal sequence
% $latex ( x_1 , \ldots , x_N )$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
% 	x_k = x\_out(:, k)
% \]$$
% and $icode x_out$$ is a two dimensional array with size $latex n \times N$$.
% 
% $head u_out$$
% The result $icode u_out$$ contains the Lagrange multiplier sequence 
% corresponding to $icode x_out$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
% 	u_k = u\_out(:, k)
% \]$$
% and $icode u_out$$ is a two dimensional array with size 
% $latex \ell \times N$$.
% The pair $icode x_out$$, $icode u_out$$ satisfy the
% $cref/first order conditions/ckbs/First Order Conditions/$$.
%
% $head info$$
% The result $icode info$$ is a matrix with each row corresponding
% to an iteration of the algorithm. 
%
% $subhead Objective$$
% The value $latex info(q, 1)$$ is the value of the objective function
% corresponding to the feasible solution found at the
% $th q$$ iteration.
%
% $subhead Step Size$$
% The value $latex info(q, 2)$$ is the step size at the $th q$$ iteration
% (as a factor of the step that solves the corresponding affine sub-problem).
%
% $subhead First Order Conditions$$
% Let $latex u_k$$ denote the dual variable vector components for
% time point $latex k$$ and at iteration $latex q$$.  
% Let $latex d_k$$ be the partial derivative of 
% $latex S ( x_1 , \ldots , x_N )$$
% with respect to $latex x_k$$ at iteration $latex q$$.
% $latex \[
% \begin{array}{rcl}
% info(q, 3) & = & \max \left|  ( B_k^T * u_k + d_k )_j        \right| 
% \\
% info(q, 4) & = & \max \left\{ (u_k)_i * (- b_k - B_k * x_k)_i \right\}
% \end{array}
% \] $$
% where the maximums are with respect to
% $latex j = 1 , \ldots , n$$,
% $latex i = 1 , \ldots , \ell$$, and
% $latex k = 1 , \ldots , N$$.
% 
% $head Example$$
% The file $cref ckbs_ok.m$$ contains an example and test of 
% $code ckbs$$.
% It returns true if $code ckbs$$ passes the test
% and false otherwise.
%
% $childtable#
%	ckbs_ok.m#
%	ckbs_blkdiag_mul.m#
%	ckbs_blkdiag_mul_t.m#
%	ckbs_sumsq_obj.m#
%	ckbs_sumsq_grad.m#
%	ckbs_sumsq_hes.m#
%	ckbs_tridiag_solve.m#
%	ckbs_newton_step.m#
%	ckbs_affine.m#
%	example/all_ok.m
% #$$
%
% $end
% ---------------------------------------------------------------------------
function [x_out, u_out, info] = ckbs(g_fun, h_fun, ...
	epsilon, x_in, z, b, db, qinv, rinv)
%
% size of problem
n     = size(x_in, 1);
N     = size(x_in, 2);
m     = size(z, 1);
ell   = size(b, 1);
%
% zero vector used in call to g_fun
x0    = zeros(n, 1);
%
% for checking sizes
[hk, Hk] = feval(h_fun, 1, x_in(:,1));
[gk, Gk] = feval(g_fun, 1, x0);
%
% check sizes
if n ~= size(db,2) | n ~= size(qinv,1) | n ~= size(qinv,2) | ...
   n ~= size(Hk,2) | n ~= size(gk,1) | n ~= size(Gk,1) | n ~= size(Gk,2)
	n
	size(db,2)
	size(qinv,1)
	size(qinv,2)
	size(Hk,2)
	size(gk,1)
	size(Gk,1)
	size(Gk,2)
	error('ckbs:  argument sizes do not agree');
end
if m ~= size(rinv,1) | m ~= size(rinv,2) | m ~= size(hk,1) | m ~= size(Hk,1)
	m
	size(rinv,1)
	size(rinv,2)
	size(hk,1)
	size(Hk,1)
	error('ckbs:  argument sizes do not agree');
end
if ell ~= size(db,1)
	ell
	size(db,1)
	error('ckbs:  argument sizes do not agree');
end
if N ~= size(z,2) | N ~= size(b,2) | N ~= size(db,3) | ...
   N ~= size(qinv,3) | N ~= size(rinv,3)
	N
	size(z,2)
	size(b,2)
	size(db,3)
	size(qinv,3)
	size(rinv,3)
	error('ckbs:  argument sizes do not agree');
end
%
% compare derivatives computed by h_fun and g_fun with central differences
check_calculation_of_derivatives = false;
if check_calculation_of_derivatives
	step = 1e-5;
	[h1, H1] = feval(h_fun, 1, x_in(:,1));
	[g1, G1] = feval(g_fun, 2, x_in(:,1));
	disp('analytic , numerical, difference');
	for j = 1 : n
		x1     = x_in(:,1);
		x1(j)  = x_in(j,1) + step;
		hp     = feval(h_fun, 1, x1);
		gp     = feval(g_fun, 2, x1);
		x1(j)  = x_in(j,1) - step;
		hm     = feval(h_fun, 1, x1);
		gm     = feval(g_fun, 2, x1);
		analytic   = [ H1(:, j) ; G1(:, j) ];
		numerical  = [ (hp - hm) ; (gp - gm) ] / (2 * step);
		difference = analytic - numerical;
		disp( [analytic , numerical , difference ] );
	end 
	keyboard('end derivative check')
end
%
% zero step corresponding to affine sub-problem
x_zero = zeros(n, N);
%
% vector version of b
b_vec = reshape(b, ell * N, 1);
%
% norm of db
if ell == 0
	db_norm = 0;
else
	db_norm = sum( sum( sum( abs(db) ) ) );
end
%
% current values evaluae h and g at x_in
x_cur  = x_in;
h_cur  = zeros(m, N);
g_cur  = zeros(n, N);
dh_cur = zeros(m, n, N);
dg_cur = zeros(n, n, N);
h      = zeros(m, N);
g      = zeros(n, N);
%
% begin the Gauss-Newton iteration
info     = zeros(0, 4);
converge = false;
itr      = 0;
while ( ~ converge ) & (itr < 100 )
	itr = itr + 1;
	%
	% value of b, g, h, db, dg, and gh at current x_cur
	xk1 = x0;
	for k = 1 : N
		xk             = x_cur(:, k);
        	[hk, Hk]       = feval(h_fun, k, xk);
                [gk, Gk]       = feval(g_fun, k, xk1);
        	dh_cur(:,:, k) = Hk;
        	dg_cur(:,:, k) = Gk;
		%
		% adjust function values so affine approximation is right
        	h_cur(:, k)    = hk - Hk * xk;
        	g_cur(:, k)    = gk - Gk * xk1;
		%
		xk1 = xk;
	end
	h_norm = sum( sum( sum( abs(h_cur) ) ) );
	g_norm = sum( sum( sum( abs(g_cur) ) ) );
	scale  = max( h_norm, g_norm );
	scale  = max( scale, db_norm);
	delta  = 1e-13 * scale;
	%
	% objective function value corresponding to x_cur
	obj_cur = ckbs_sumsq_obj(x_cur, z, ...
		g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
	%
	[x_new, u_new, affine_info] = ckbs_affine(delta, x_in, z, ...
		b, g_cur, h_cur, db, dg_cur, dh_cur, qinv, rinv);
	%
        % affine approximate objective function at x_cur + x_step
        obj_affine = ckbs_sumsq_obj(x_new, z, ...
                g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
	diff_affine = obj_cur - obj_affine;
	%
	% check for convergence
	d_cur  = ckbs_sumsq_grad(x_cur, z, ...
		g_cur, h_cur, dg_cur, dh_cur, qinv, rinv);
	d_vec = reshape(d_cur, n * N, 1);
	if ell > 0
		u_vec = reshape(u_new, ell * N, 1);
		x_vec = reshape(x_cur, n * N, 1);
		s_vec = - b_vec - ckbs_blkdiag_mul(db, x_vec);
		%
        	% check for feasibility
        	if min(s_vec) <= 0 | min(u_vec) <= 0
                	error('ckbs: current iterate is not feasible');
        	end
		Bt_u_d   = ckbs_blkdiag_mul_t(db, u_vec) + d_vec;
		e1       = max( abs( Bt_u_d ) );
		e2       = max( u_vec .* s_vec );
	else
		e1       = max( abs(d_vec) );
		e2       = 0;
	end
	converge =  ( e1 < epsilon) & ( e2 < epsilon );
	%
	% line search
	if converge
		lambda = 0;
	else
		lambda = 2;
		if diff_affine < 0
			error('ckbs: affine objective increased');
		end
		x_step    = x_new - x_cur;
	end
	ok        = converge;
	kount     = 0;
	max_kount = 10;
	while (~ ok ) & (kount < max_kount)
		kount  = kount + 1;
		lambda = lambda / 2;
		%
		% argument corresponding to this step
		x = x_cur + lambda .* x_step;
		%
		% value of g, corresponding to x
		xk1 = x0;
		for k = 1 : N
			xk       = x(:, k);
        		hk       = feval(h_fun, k, xk);
                	gk       = feval(g_fun, k, xk1);
        		Hk       = dh_cur(:,:, k);
        		Gk       = dg_cur(:,:, k);
        		h(:, k)  = hk - Hk * xk;
        		g(:, k)  = gk - Gk * xk1;
			xk1 = xk;
		end
		%
		% true objective corresponding to xcur + lambda * step
		obj_true = ckbs_sumsq_obj(x, z, ...
			g, h, dg_cur, dh_cur, qinv, rinv);
		%
		diff_true = obj_cur - obj_true;
		ok  = diff_true >= lambda * diff_affine / 2;
	end
	if ~ok 
		error('ckbs: line search failed');
		return;
	end
	% values corresponding to this iteration
	info_itr = [ obj_cur, lambda , e1 , e2 ]
	info     = [ info ; info_itr ];
	%
	% value for next iteration
	x_cur   = x;
end
x_out = x_cur;
u_out = u_new;
if ~ converge
        error('ckbs: did not converge');
end
return
end
