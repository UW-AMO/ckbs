% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Aravkin:    saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_t_obj$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_t_obj
%       qinv
%       rinv
%       dh
%       dg
%       params
%       pos
%       vel
%       df
%       inds
%       res
%       gk
%       feval
%       dt
% $$
%
% $section Student's t Sum of Squares Objective$$
%
% $index ckbs_t_obj$$
% $index t_obj$$
%
% $index objective, Student's t$$
% $index Student's t, objective$$
%
% $head Syntax$$
% $codei/[/obj/] = ckbs_t_obj(/
%       x/, /z/, /g_fun/, /h_fun/, /qinv/, /rinv/, /params/)/$$
%
% $head Purpose$$
% This routine computes the value of the
% Student's t objective function.
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
% \\
% \end{array}
% \] $$
% where the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive definite and
% $latex x_0$$ is the constant zero.
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
% The argument $icode g_fun$$ is a function handle to the process model. 
%
% $head h_fun$$
% The argument $icode h_fun$$ is a function handle to the measurement
% model.
%
% $head qinv$$
% The argument $icode qinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       Q_k^{-1} = qinv(:,:, k)
% \]$$
% and $icode qinv$$ has size $latex n \times n \times N$$.
%
% $head rinv$$
% The argument $icode rinv$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       R_k^{-1} = rinv(:,:, k)
% \]$$
% and $icode rinv$$ has size $latex m \times m \times N$$.
%
% $head obj$$
% The result $icode obj$$ is a scalar given by
% $latex \[
%       obj = S ( x_1 , \ldots , x_N )
% \] $$
%
% $children#
%       example/t_obj_ok.m
% #$$
%
% $head Example$$
% The file $cref t_obj_ok.m$$ contains an example and test of
% $code ckbs_t_obj$$.
% It returns true if $code ckbs_t_obj$$ passes the test
% and false otherwise.
%
% $end
% -------------------------------------------------------------------------
% ---

function [obj] = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params)


m    = size(z, 1);
n    = size(x, 1);
N    = size(x, 2);

res_z = zeros(m, N);
for k = 1 : N
    xk    = x(:, k);
    res_z(:,k)    = z(:,k) - feval(h_fun, k, xk);
end


res_x = zeros(n, N);
xk = zeros(n, 1);
for k = 1 : N
	xkm   = xk;
	xk    = x(:, k);
	gk    = feval(g_fun, k, xkm);
	res_x(:,k)  = xk - gk;
end


%----------------------------------------------------------------------
% must be included in params
df_proc = params.df_proc; 
df_meas = params.df_meas;

inds_proc_st = params.inds_proc_st;        % which indices to apply student's t to
inds_proc_sq = setdiff(1:n, inds_proc_st); 

inds_meas_st = params.inds_meas_st; % which measurement indices to apply student's t to
inds_meas_sq = setdiff(1:m, inds_meas_st);
%----------------------------------------------------------------------



if(~isempty(inds_proc_st))
   res_x_st = res_x(inds_proc_st, :);   
   qinv_st = qinv(inds_proc_st, inds_proc_st, :);
else
   res_x_st = zeros(1,N); 
   qinv_st = zeros(1,1,N);
end

if(~isempty(inds_proc_sq))
   res_x_sq = res_x(inds_proc_sq, :);   
   qinv_sq = qinv(inds_proc_sq, inds_proc_sq, :);
else
   res_x_sq = zeros(1,N); 
   qinv_sq = zeros(1,1,N);
end
if(~isempty(inds_meas_st))
   res_z_st = res_z(inds_meas_st, :);  
   rinv_st = rinv(inds_meas_st, inds_meas_st,:);
else
   res_z_st = zeros(1,N);
   rinv_st = zeros(1,1,N);
end

if(~isempty(inds_meas_sq))
   res_z_sq = res_z(inds_meas_sq, :);
   rinv_sq = rinv(inds_meas_sq, inds_meas_sq,:);
else
   res_z_sq = zeros(1,N); 
   rinv_sq = zeros(1,1,N);
end


obj  = 0;

for k = 1 : N
    
    
    qinv_st_k = qinv_st(:, :, k);  
    qinv_sq_k = qinv_sq(:, :, k);
    rinv_st_k = rinv_st(:, :, k);
    rinv_sq_k = rinv_sq(:, :, k);

    res_x_st_k = res_x_st(:,k);
    res_x_sq_k = res_x_sq(:,k);
    res_z_st_k = res_z_st(:,k);
    res_z_sq_k = res_z_sq(:,k);
      

    obj   = obj + 0.5*(df_meas) * log(1 + (1/df_meas)*res_z_st_k' * rinv_st_k * res_z_st_k) ...
                + 0.5*(df_proc) * log(1 + (1/df_proc)*res_x_st_k' * qinv_st_k * res_x_st_k) ...
                + 0.5 * res_z_sq_k' * rinv_sq_k * res_z_sq_k ...
                + 0.5 * res_x_sq_k' * qinv_sq_k * res_x_sq_k;

end




return
