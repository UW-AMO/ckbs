% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_t_hess$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_t_hess
%       dg
%       dh
%       qinv
%       rinv
%       params
% $$
%
% $section Student's t Hessian$$
%
% $index ckbs_t_hess$$
% $index t_hess$$
%
% $index hessian, of Student's t objective$$
% $index Students t objective, hessian$$
%
% $head Syntax$$
% $codei/[/grad/] = ckbs_t_hess(/
%       x/, /z/, /g_fun/, /h_fun/, /qinv/, /rinv/, /params/)/$$
%
% $head Purpose$$
% This computes the gradient of the
% of the general Student's t objective.
%
% $head Notation$$
% The affine Kalman-Bucy smoother residual sum of squares is defined by
% $latex \[
% \begin{array}{rcl}
% S ( x_1 , \ldots , x_N ) & = & \sum_{k=1}^N S_k ( x_k , x_{k-1} ) \\
% S_k ( x_k , x_{k-1} )    & = &
% \frac{1}{2}
% ( z_k - h_k - H_k * x_k )^\R{T} * R_k^{-1} * ( z_k - h_k - H_k * x_k )
% \\
% & + &
% \frac{1}{2}
% ( x_k - g_k - G_k * x_{k-1} )^\R{T} * Q_k^{-1} * ( x_k - g_k - G_k * x_{k-1} )
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
%
% $head Gradient$$
% We define $latex Q_{N+1}$$ to be the $latex n \times n$$ identity
% matrix and $latex G_{N+1}$$ to be zero,
% $latex \[
% \begin{array}{rcl}
% \nabla_k S_k^{(1)} ( x_k , x_{k-1} )
% & = &  H_k^\R{T} * R_k^{-1} * ( h_k + H_k * x_k - z_k )
%   +    Q_k^{-1} * ( x_k - g_k - G_k * x_{k-1} )
% \\
% \nabla_k S_{k+1}^{(1)} ( x_{k+1} , x_k )
% & = & G_{k+1}^\R{T} * Q_{k+1}^{-1} * ( g_{k+1} + G_{k+1} * x_k  - x_{k+1} )
% \end{array}
% \] $$
% It follows that the gradient of the
% affine Kalman-Bucy smoother residual sum of squares is
% $latex \[
% \begin{array}{rcl}
% \nabla S ( x_1 , \ldots , x_N )
% & = &
% \left( \begin{array}{c}
%       d_1 \\ \vdots \\ d_N
% \end{array} \right)
% \\
% d_k & = & \nabla_k S_k^{(1)}     ( x_k , x_{k-1} )
%       +   \nabla_k S_{k+1}^{(1)} ( x_{k+1} , x_k )
% \end{array}
% \] $$
% where $latex S_{N+1} ( x_{N+1} , x_N )$$ is defined as
% identically zero.
%
% $head x$$
% The argument $icode x$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       x_k = x(:, k)
% \]$$
% and $icode x$$ has size $latex n \times N$$.
%
% $head z$$
% The argument $icode z$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       z_k = z(:, k)
% \]$$
% and $icode z$$ has size $latex m \times N$$.
%
% $head g_fun$$
% The argument $icode g_fun$$ is a function handle for the
% process model.
%
% $head h_fun$$
% The argument $icode h_fun$$ is a function handle for the 
% measurement model. 
%
% $head qinv$$
% The argument $icode qinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       Q_k^{-1} = qinv(:,:,k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
%
% $head rinv$$
% The argument $icode rinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       R_k^{-1} = rinv(:,:,k)
% \]$$
% and $icode rinv$$ has size $latex m \times m \times N$$.
%
% $head params$$
% The argument $icode params$$ is a structure containing the 
% requisite parameters. 
%
% $head D$$
% The result $icode D$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       D_k = hess(:, :, k)
% \]$$
% and $icode D$$ has size $latex n\times n \times N$$.
%
% $head A$$
% The result $icode A$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       A_k = hess(:,:,k)
% \]$$
% and $icode A$$ has size $latex n\times n \times N$$.

% $children#
%       example/t_hess_ok.m
% #$$
%
% $head Example$$
% The file $cref t_hess_ok.m$$ contains an example and test of
% $code ckbs_t_hess$$.
% It returns true if $code ckbs_t_hess$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------
function [D A] = ckbs_t_hess(x, z, g_fun, h_fun, qinv, rinv, params)


%
% sizes for this problem
m = size(z, 1);
n = size(x, 1); 
N = size(x, 2);
%
% initialize return values
% dimension return values
D  = zeros(n, n, N);
A  = zeros(n, n, N);

Gx = zeros(n, N);
xk = zeros(n, 1);
for k = 1 : N
	xkm   = xk;
	xk    = x(:, k);
	gk    = feval(g_fun, k, xkm);
	Gx(:,k)  = xk - gk;
end

Hx = zeros(m, N);
for k = 1 : N
    xk    = x(:, k);
    Hx(:,k)    = feval(h_fun, k, xk);
end

%

%----------------------------------------------------------------------
% must be included in params
df_proc = params.df_proc; 
df_meas = params.df_meas;

inds_proc_st = params.inds_proc_st;        % which indices to apply student's t to
inds_proc_sq = setdiff(1:n, inds_proc_st); 

inds_meas_st = params.inds_meas_st; % which measurement indices to apply student's t to
inds_meas_sq = setdiff(1:m, inds_meas_st);
%----------------------------------------------------------------------

Gk    = zeros(n, n);
Qkinv = eye(n, n);

if(~isempty(inds_proc_st))
   Gk_st = Gk(inds_proc_st, inds_proc_st);
   Qkinv_st = Qkinv(inds_proc_st, inds_proc_st);
else
   Gk_st = 0;
   Qkinv_st = 0;
end

if(~isempty(inds_proc_sq))
   Gk_sq = Gk(inds_proc_sq, inds_proc_sq);
   Qkinv_sq = Qkinv(inds_proc_sq, inds_proc_sq);
else
   Gk_sq = 0; 
   Qkinv_sq = 0;
end

for k = N : -1 : 1
	Gk1_st     = Gk_st;
    Gk1_sq     = Gk_sq;
	Qk1inv_st  = Qkinv_st;
    Qk1inv_sq  = Qkinv_sq;
	%
	xk      = x(:, k);
    [gk, Gk]           = feval(g_fun, k, xk);
	[hk, Hk]            = feval(h_fun, k, xk);   
	Qkinv   = qinv(:,:, k);
	Rkinv   = rinv(:,:, k);
    if( k == N)
       xk1res = zeros(n,1); 
    else
       xk1res = Gx(:, k+1);
    end 
    xkres = Gx(:, k);
    zkres = Hx(:, k) - z(:,k);

    
if(~isempty(inds_proc_st))
   xk1res_st = xk1res(inds_proc_st, :);
   xkres_st = xkres(inds_proc_st, :);
   Gk_st = Gk(inds_proc_st, inds_proc_st);
   Qkinv_st = Qkinv(inds_proc_st, inds_proc_st, :);
else
   xk1res_st = 0;  
   xkres_st = 0; 
   Gk_st = 0;
   Qkinv_st = 0;
end

if(~isempty(inds_proc_sq))
   Gk_sq    = Gk(inds_proc_sq, inds_proc_sq);
   Qkinv_sq = Qkinv(inds_proc_sq, inds_proc_sq, :);
else
   Gk_sq = 0;
   Qkinv_sq = 0;
end
if(~isempty(inds_meas_st))
   zkres_st = zkres(inds_meas_st, :);
   Hk_st    = Hk(inds_meas_st, :);
   Rkinv_st = Rkinv(inds_meas_st, inds_meas_st,:);
else
   zkres_st = 0;
   Hk_st = 0;
   Rkinv_st = 0;
end

if(~isempty(inds_meas_sq))
   Hk_sq = Hk(inds_meas_sq, :);
   Rkinv_sq = Rkinv(inds_meas_sq, inds_meas_sq,:);
else
   Hk_sq = 0;
   Rkinv_sq = 0;
end
    
    % Compute weighted factors
    sk      = (df_meas)/(df_meas + zkres_st'  * Rkinv_st  * zkres_st);
    rk      = (df_proc)/(df_proc + xkres_st'  * Qkinv_st  * xkres_st);
    rk1     = (df_proc)/(df_proc + xk1res_st' * Qk1inv_st * xk1res_st);
 
    
    if(~isempty(inds_proc_st))
        D(inds_proc_st,inds_proc_st, k) = rk * Qkinv_st ...
             + rk1 * Gk1_st' * Qk1inv_st * Gk1_st;
        A(inds_proc_st, inds_proc_st, k) = -rk*(Qkinv_st * Gk_st);
    end
    if(~isempty(inds_proc_sq))
        D(inds_proc_sq,inds_proc_sq, k) = Qkinv_sq ...
             + Gk1_sq' * Qk1inv_sq * Gk1_sq;
        A(inds_proc_sq, inds_proc_sq, k) = -Qkinv_sq * Gk_sq;
    end


    D(:,:,k) = D(:,:,k) + sk * Hk_st' * Rkinv_st * Hk_st ...
        + Hk_sq' * Rkinv_sq * Hk_sq;
    
    
    

%     D(:,:,k) = sk * Hk_st' * Rkinv_st * Hk_st ...
%              + rk1 * Qkinv_st ...
%              + rk * Gk1_st' * Qk1inv_st * Gk1_st...
%              + Hk_sq' * Rkinv_sq * Hk_sq...
%              + Qkinv_sq ...
%              + Gk1_sq' * Qk1inv_sq * Gk1_sq;
% 	A(:,:,k) = -rk*(Qkinv_st * Gk_st) - Qkinv_sq * Gk_sq;

end
return
