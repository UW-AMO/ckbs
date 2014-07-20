% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_t_general$$ $newlinech %$$
% $spell
%       itr
%       obj
%       ckbs
%       qinv
%       rinv
%       feval
%       fk
%       gk
%       hk
%       xk
%       fu
%       aff
%       params
%       df
%       dif
%       inds
%       indicies
% $$
%
% $section The General Student's t Smoother$$
%
% $index ckbs_t_general$$
% $index Student's t, Kalman-Bucy Smoother$$
% $index smoother, Student's t$$
% $index Kalman, Student's t smoother$$
%
% $head Syntax$$
% $codei/[/x_out/, /info/] = ckbs_t_general(/
%       /g_fun/, /h_fun/, ...
%       /max_itr/, /epsilon/, /x_in/, /z/, /qinv/, /rinv/, /params/)/$$
%
% $head Purpose$$
% This routine minimizes the objective that results when 
% Student's t errors are used for user specified components of the 
% process and measurement model. 
% Process and measurement models must be passed in as function handles, 
% so even in the affine case function handles, rather than vectors and 
% matrices, must be provided. 
%
% 
%
% $head Notation$$
% The general objective solved by the smoother is given by 
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = &
% \rho_k^p\left([ z_k - h_k ( x_k ) ]^\R{T} * R_k^{-1} * [ z_k - h_k ( x_k
% ) ]\right)
% \\
% & + &
% \rho_k^m\left([ x_k - g_k ( x_{k-1} ) ]^\R{T} * Q_k^{-1} * [ x_k - g_k (
% x_{k-1} ) ]\right)
% \\
% \end{array}
% \] $$
% where $latex \rho_k^p$$ and $latex \rho_k^m$$ are either quadratic, Student's t, 
% or component-wise mixtures of these functions, 
% depending on user-specified indices.  The matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
% Note that the process model $latex g_fun$$ should provide an initial state estimate for $latex k=1$$, 
% and $latex Q_1$$ is the corresponding covariance.
%
% $head Problem$$
% The full general Student's t problem is given by 
% $latex \[
% \begin{array}{rll}
% {\rm minimize}
%       & S( x_1 , \ldots , x_N )
%       & {\rm w.r.t.} \; x_1 \in \B{R}^n , \ldots , x_N \in \B{R}^n
% \end{array}
% \] $$
%
%
%
% $head g_fun$$
% The $code ckbs_nonlinear$$ argument $icode g_fun$$
% is a function that supports both of
% the following syntaxes
% $codei/
%       [/gk/] = feval(/g_fun/, /k/, /xk1/)
%       [/gk/, /Gk/] = feval(/g_fun/, /k/, /xk1/)
% /$$
%
% $subhead k$$
% The $icode g_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$.
% The case $latex k = 1$$ serves to specify the initial state estimate.
%
% $subhead xk1$$
% The $icode g_fun$$ argument $icode xk1$$ is an column vector with
% length $latex n$$. It specifies the state vector at time index $icode k$$
% $latex \[
%       xk1 = x_{k-1}
% \] $$.
% In the case $latex k = 1$$, the value of $icode xk1$$ does not matter.
%
% $subhead gk$$
% The $icode g_fun$$ result $icode gk$$ is an column vector with
% length $latex n$$ and
% $latex \[
%       gk = g_k ( x_{k-1} )
% \] $$
% In the case $latex k = 1$$, the value of $icode gk$$ is the initial state
% estimate at time index $icode k$$.
%
% $subhead Gk$$
% If the $icode g_fun$$ result $icode Gk$$ is present in the syntax,
% it is the $latex n \times n$$ matrix given by
% and
% $latex \[
%       Gk = \partial_{k-1} g_k ( x_{k-1} )
% \] $$
% In the case $latex k = 1$$, the value of $icode Gk$$ is the zero matrix;
% i.e., $latex g_k$$ does not depend on the value of $latex x_{k-1}$$.
%
% $head h_fun$$
% The $code ckbs_nonlinear$$ argument $icode h_fun$$
% is a function that supports both of the
% following syntaxes
% $codei/
%       [/hk/] = feval(/h_fun/, /k/, /xk/)
%       [/hk/, /Hk/] = feval(/h_fun/, /k/, /xk/)
% /$$
%
% $subhead k$$
% The $icode h_fun$$ argument $icode k$$ is an integer with
% $latex 1 \leq k \leq N$$.
%
% $subhead xk$$
% The $icode h_fun$$ argument $icode xk$$ is an column vector with
% length $latex n$$. It specifies the state vector at time index $icode k$$
% $latex \[
%       xk = x_k
% \] $$.
%
% $subhead hk$$
% The $icode h_fun$$ result $icode hk$$ is an column vector with
% length $latex m$$ and
% $latex \[
%       hk = h_k ( x_k )
% \] $$
%
% $subhead Hk$$
% If the $icode h_fun$$ result $icode Hk$$ is present in the syntax,
% it is the $latex m \times n$$ matrix given by
% and
% $latex \[
%       Hk = \partial_k h_k ( x_k )
% \] $$
%
% $head max_itr$$
% The integer scalar $icode max_itr$$ specifies the maximum number of
% iterations of the algorithm to execute. It must be greater than or
% equal to zero. Note that if it is zero, the first row of the
% $cref/info/ckbs_nonlinear/info/$$ return value will still be computed.
% This can be useful for deciding what is a good value for the argument
% $cref/epsilon/ckbs_nonlinear/epsilon/$$.
%
% $head epsilon$$
% The $code ckbs_nonlinear$$ argument $icode epsilon$$ is a positive scalar.
% It specifies the convergence
% criteria value; i.e.,
% $latex \[
%       \varepsilon = epsilon
% \] $$
%
% $head x_in$$
% The $code ckbs_nonlinear$$ argument $icode x_in$$
% is a two dimensional array with size $latex n \times N$$.
% It specifies a sequence of state values; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x\_in (:, k) = x_k
% \] $$
% The closer the initial state sequence is to the solution
% the faster, and more likely, the $code ckbs_nonlinear$$ will converge.
% The initial state sequence need not be feasible; i.e.
% it is not necessary that
% $latex \[
%       f_k ( x_k ) \leq 0
% \] $$
% for all $latex k$$.
%
% $head z$$
% The $code ckbs_nonlinear$$ argument $icode z$$ is a two dimensional array
% with size $latex m \times N$$.
% It specifies the sequence of measurement vectors; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       z(:, k) = z_k
% \]$$
%
% $head qinv$$
% The $code ckbs_nonlinear$$ argument $icode qinv$$
% is a three dimensional array
% with size $latex n \times n \times N$$.
% It specifies the inverse of the variance of the measurement noise; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       qinv(:,:, k) = Q_k^{-1}
% \]$$
% In the case $latex k = 1$$, the value of $latex Q_k$$ is the variance
% of the initial state estimate (see $cref/g_fun/ckbs_nonlinear/g_fun/$$.
%
% $head rinv$$
% The $code ckbs_nonlinear$$ argument $icode rinv$$
% is a three dimensional array,
% with size $latex m \times m \times N$$.
% it specifies the inverse of the variance of the transition noise; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       rinv(:,:, k) = R_k^{-1}
% \]$$
% It is ok to signify a missing data value by setting the corresponding
% row and column of $icode rinv$$ to zero. This also enables use
% to interpolate the state vector $latex x_k$$ to points where
% there are no measurements.
%
% $head params$$
% The $icode params$$ structure specifies auxiliary parameters required by 
% the smoother. The structure is passed to the objective, gradient, and
% hessian functions; some parameter values are required for these to work
% correctly (see below). 
%
% $subhead level$$
% The $icode level$$ argument sets the output level. If $icode level$$
% is greater than 1, numerical derivatives are computed. 
%
% $subhead df_proc$$
% The $icode dif_proc$$ is the Student's t degree of freedom parameter used
% for any Student's t process component. 
%
% $subhead df_meas$$
% The $icode dif_meas$$ is the Student's t degree of freedom parameter used
% for any Student's t measurement component. 
%
% $subhead inds_proc_st$$
% The $icode inds_proc_st$$ lists indicies of the state (process residual)
% that are to be modeled using Student's t. Setting this structure to empty 
% $latex []$$ means quadratic model is used for entire state. 
%
% $subhead inds_meas_st$$
% The $icode inds_meas_st$$ lists indicies of the measurement residual
% that are to be modeled using Student's t. Setting this structure to empty 
% $latex []$$ means quadratic model is used for entire measurement. 
%
% $head x_out$$
% The result $icode x_out$$
% is a two dimensional array with size $latex n \times N$$.
% $latex ( x_1 , \ldots , x_N )$$.
% It contains an approximation for the optimal sequence; i.e.,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x\_out(:, k) = x_k
% \]$$
%
%
% $head info$$
% The result $icode info$$ is a matrix with each row corresponding
% to an iteration of the algorithm.
% Note that $code ckbs_nonlinear$$ has satisfied the convergence condition if
% and only if
% $codei%
%       all( %info%(end, 1) <= %epsilon% )
% %$$
%or 
% $codei%
%       all( %info%(end, 3) <= %epsilon% )
% %$$
%
% $children#
%       example/t_general_ok.m#
%       example/t_general_noisy_jump.m
% #$$
%
% $head Example$$
%
% $subhead Simple$$
% The file $cref t_general_ok.m$$ contains a simple example
% and test of $code ckbs_t_general$$.
% It returns true if $code ckbs_t_general$$ passes the test 
% (meaning algorithm converges to stationary point)
% and false otherwise.
%
% $end

function [xOut, info] = ckbs_t_general(g_fun, h_fun, ...
	max_itr, epsilon, x_in, z, qinv, rinv, params)

level = params.level; 

if nargin ~= 9
	error('ckbs_t_general: improper number of input arguments');
end
n         = size(x_in, 1);

[gk, Gk] = feval(g_fun, 1, zeros(n, 1));
[hk, Hk] = feval(h_fun, 1, x_in(:,1));
%
% get other sizes of problem
N     = size(x_in, 2);
m     = size(z, 1);
%
% check n in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
if n ~= size(x_in,1) || n ~= size(qinv,1) || n ~= size(qinv,2) || ...
   n ~= size(gk,1) || n ~= size(Gk,1) || n ~= size(Gk,2) || ...
   n ~= size(Hk,2)
	size(x_in,1)
	size(qinv,1)
	size(qinv,2)
	size(gk,1)
	size(Gk,1)
	size(Gk,2)
	size(Hk,2)
	error('ckbs_nonlinear: argument sizes with value n do not agree');
end
% check m in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
if m ~= size(z,1) || m ~= size(rinv,1) || m ~= size(rinv,2) || ...
   m ~= size(hk,1) || m ~= size(Hk,1)
	size(z,1)
	size(rinv,1)
	size(rinv,2)
	size(hk,1)
	size(Hk,1)
	error('ckbs_nonlinear:  argument sizes with value m do not agree');
end
% check N in following order: x_in, z, qinv, rinv, fk, Fk, gk, Gk, hk, Hk
if N ~= size(x_in,2) || N ~= size(z,2) || N ~= size(qinv,3) || N ~= size(rinv,3)
	size(x_in,2)
	size(z,2)
	size(qinv,3)
	size(rinv,3)
	error('ckbs_nonlinear:  argument sizes with value N do not agree');
end
%
% compare derivatives computed by f_fun h_fun and g_fun 
% with central difference approximations
if level >= 2
	step = 1e-5;
	[g1, G1] = feval(g_fun, 2, x_in(:,1));
	[h1, H1] = feval(h_fun, 1, x_in(:,1));
	disp('analytic , numerical, difference');
	for j = 1 : n
		x1     = x_in(:,1);
		x1(j)  = x_in(j,1) + step;
		gp     = feval(g_fun, 2, x1);
		hp     = feval(h_fun, 1, x1);
		%
		x1(j)  = x_in(j,1) - step;
		hm     = feval(h_fun, 1, x1);
		gm     = feval(g_fun, 2, x1);
		%
		analytic   = [ G1(:, j) ; H1(:, j) ];
		numerical  = [ (gp-gm) ; (hp-hm) ] / (2*step);
		difference = analytic - numerical;
		disp( [analytic , numerical , difference ] );
	end 
	keyboard('ckbs_t_general: end derivative check')
end
%
% a useful constant
zero_n     = zeros(n, 1);
%
% dimension of array
%
g_cur  = zeros(n, N);
h_cur  = zeros(m,  N);
dg_cur = zeros(n, n, N);
dh_cur = zeros(m, n, N);
x_cur  = zeros(n, N);
%
%  initialize some values
 info       = zeros(0, 4); % return convergence tracking information
% lambda     = 0;           % linear search step size
% max_affine = 50;          % maximum number of iterations in affine sub-problem
% alpha      = 0;           % initial exact penalty function parameter

obs        = n   * N;
itr        = 0;           % iteration counter


%
% initialize current value
xk1      = zero_n;
%dist_cur = 0;
for k  = 1 : N
	xk                           = x_in(:, k);
	[gk, Gk]                     = feval(g_fun, k, xk1);
	[hk, Hk]                     = feval(h_fun, k, xk);
	%
%	dist_cur                     = dist_cur + sum( max(fk, 0 ) );
	%
	% adjust so affine approximations are relative to zero (not x_cur)
	g_cur(:, k)                  = gk - Gk * xk1;
	h_cur(:, k)                  = hk - Hk * xk;
	%
	x_cur(:, k)                  = xk;
	dg_cur(:,:, k)               = Gk;
	dh_cur(:,:, k)               = Hk;
	%
	xk1                          = xk;
end

converge = 0; 

while ( ~ converge ) && (itr < max_itr)

    
    Curr_Obj = ckbs_t_obj(x_cur, z, g_fun, h_fun, qinv, rinv, params) ;   
    
	itr = itr + 1;
	%
	
    [D, A] = ckbs_t_hess(x_cur, z, g_fun, h_fun, qinv, rinv, params);

    
    % Compute gradient according to t paradigm

    %d = t_nonlin_grad(x_cur, z, g_fun, h_fun,  qinv, rinv, tdf, tdim);
    d = ckbs_t_grad(x_cur, z, g_fun, h_fun, qinv, rinv, params);
    
    dVec   = reshape(d, obs, 1);
    
    [y_new] = ckbs_tridiag_solve_b(D, A, - dVec);
    y_new = reshape(y_new, n, N);

    
   c = 0.001;
   %c = 0;
   gamma = 0.5;
   lambda_s = 2.;
   
   
   done = false;
   search_itr = 0;
   max_search_itr = 100;
   
   % WARNING! Changed y_new - y_cur to y_new
   dir = reshape(y_new, obs, 1);    
  % 'gradient norm'
  % norm(dVec,2)
  % 'normalized dirDer'
   dirDer = dVec'*dir;
  
   info(itr, 1) = dirDer;
   
   if(dirDer >=0)
    %xOut = y   ;
    %y_new = -d;
    error('positive dirDer in ckbs_t_general');
   end
   
   % don't do line search if done. 
   if(abs(dirDer) < 1e-8)
      xOut     = x_cur;
%      'ended because of small dirDer'
      return;
   end
   
   
   obj_cur = ckbs_t_obj(x_cur, z, g_fun, h_fun, qinv, rinv, params);
   info(itr, 2) = obj_cur;
   
   while(~done)&&(search_itr < max_search_itr)
       search_itr = search_itr + 1;
       lambda_s = lambda_s * gamma;
              
       % WARNING! Changed y_new - y_cur to y_new
       x_lambda_s = x_cur + lambda_s*(y_new); 
       obj_lambda_s = ckbs_t_obj(x_lambda_s, z, g_fun, h_fun, qinv, rinv, params);
       done = ((obj_lambda_s - obj_cur) <= c*lambda_s*dirDer); 
       %search_itr
       %obj_lambda_s - obj_cur
       %c*lambda_s*dVec'*dir
       
   end
   
   if(search_itr == max_search_itr)
      error('line search did not converge'); 
   end	
    
	x_cur    = x_lambda_s;
	
    
    
    % new value for objective function
	xOut     = x_cur;

    infNorm = max(dVec);
    info(itr, 3) = infNorm;
    converge = (infNorm <= epsilon);
    
    End_Obj = ckbs_t_obj(x_cur, z, g_fun, h_fun, qinv, rinv, params);
    
    d_obj = Curr_Obj - End_Obj;
    info(itr, 4) = d_obj;
    if(d_obj < epsilon*1e-4)
      %  'ended because of small d_obj'
       return 
    end
    
    
    
end
return
end
