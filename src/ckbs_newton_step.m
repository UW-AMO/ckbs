% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_newton_step$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       Hdia
%       Hlow
%       Bdia
%       tri
%       ds
%       dy
%       du
%       mu
%       Kuhn
%   tridiagonal
% $$
%
% $section Affine Constrained Kalman Bucy Smoother Newton Step$$
%
% $index ckbs_newton_step$$
% $index newton_step$$
%
% $index step, Newton$$
%
% $head Syntax$$
% $codei/[/ds/, /dy/, /du/] = ckbs_newton_step(/
%       mu/, /s/, /y/, /u/, /b/, /d/, /Bdia/, /Hdia/, /Hlow/)/$$
%
% $head Purpose$$
% This routine computes one step of Newton's method for solving
% the non-linear Kuhn-Tucker conditions for
% the $latex \mu$$-relaxed affine constrained Kalman-Bucy smoother problem.
%
%
% $head Problem$$
% Given
% $latex \mu \in \B{R}_+$$,
% $latex H \in \B{R}^{p \times p}$$,
% $latex d \in \B{R}^p$$,
% $latex b \in \B{R}^r$$, and
% $latex B \in \B{R}^{r \times p}$$,
% the $latex \mu$$-relaxed affine constrained Kalman-Bucy smoother problem is:
% $latex \[
% \begin{array}{rl}
% {\rm minimize} & \frac{1}{2} y^\R{T} H y + d^\R{T} y
% - \mu \sum_{i=1}^r \log(s_i)
% \; {\rm w.r.t} \; y \in \B{R}^p \; , \; s \in \B{R}_+^r
% \\
% {\rm subject \; to} & s + b + B y  = 0
% \end{array}
% \] $$
% In addition, $latex H$$ is symmetric block tri-diagonal with each block of
% size $latex n \times n$$ and
% $latex B$$ is block diagonal with each block of size $latex m \times n$$
% (there is an integer $latex N$$ such that
% $latex p = n * N$$ and $latex r = m * N$$).
%
%
% $subhead Lagrangian$$
% We use $latex u \in \B{R}^r$$
% to denote the Lagrange multipliers corresponding to the constraint equation.
% The corresponding Lagrangian is
% $latex \[
% L(y, s, u)  =
% \frac{1}{2} y^\R{T} H y + d^\R{T} y
% - \mu \sum_{i=1}^r \log(s_i)
% + u^\R{T} (s + b + B y)
% \] $$
% The partial gradients of the Lagrangian are given by
% $latex \[
% \begin{array}{rcl}
% \nabla_y L(y, s, u ) & = & H y + B^\R{T} u + d  \\
% \nabla_s L(y, s, u ) & = & u - \mu / s \\
% \nabla_u L(y, s, u ) & = & s + b + B y \\
% \end{array}
% \] $$
% where $latex \mu / s $$ is the component by component division of
% $latex \mu $$ by the components of the $latex s$$.
% Note, from the second equation, that we only need consider
% $latex u \geq 0$$ because $latex s \geq 0$$.
%
% $subhead Kuhn-Tucker Conditions$$
% We use $latex D(s)$$
% to denote the diagonal matrix with $latex s$$
% along its diagonal and $latex 1_r$$
% to denote the vector, of length $latex r$$ with all its components
% equal to one.
% We define
% $latex F : \B{R}^{r + p + r} \rightarrow \B{R}^{r + p + r}$$
% by
% $latex \[
% F(s, u, y)
% =
% \left(
% \begin{array}{c}
% s + b + B y       \\
% D(s) D(u) 1_r - \mu 1_r\\
% H y + B^\R{T} u + d
% \end{array}
% \right)
% \] $$
% The Kuhn-Tucker conditions for a solution of the
% $latex \mu$$-relaxed constrained affine Kalman-Bucy smoother problem is
% $latex F(s, u, y) = 0 $$.
%
% $subhead Newton Step$$
% Given a value for $latex (s, u, y)$$, the Newton step
% $latex ( \Delta s^\R{T} , \Delta u^\R{T} , \Delta y^\R{T} )^\R{T}$$ solves the problem:
% $latex \[
% F^{(1)} (s, u, y)
% \left( \begin{array}{c} \Delta s \\ \Delta u \\ \Delta y \end{array} \right)
% =
% - F(s, u, y)
% \] $$
%
% $head mu$$
% The argument $icode mu$$ is a positive scalar specifying the
% relaxation parameter $latex \mu$$.
%
% $head s$$
% The argument $icode s$$ is a column vector of length $latex r$$.
% All the elements of $icode s$$ are greater than zero.
%
% $head y$$
% The argument $icode y$$ is a column vector of length $latex p$$
%
% $head u$$
% The argument $icode u$$ is a column vector of length $latex r$$.
% All the elements of $icode s$$ are greater than zero.
%
% $head b$$
% The argument $icode b$$ is a column vector of length $latex r$$.
%
% $head d$$
% The argument $icode d$$ is a column vector of length $latex p$$
%
% $head Bdia$$
% The argument $icode Bdia$$ is an $latex m \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex B_k \in \B{R}^{m \times n}$$ by
% $latex \[
%       B_k = Bdia(:, :, k)
% \] $$
%
% $head B$$
% The matrix $latex B$$ is defined by
% $latex \[
% B
% =
% \left( \begin{array}{cccc}
% B_1 & 0      & 0      &           \\
% 0   & B_2    & 0      & 0         \\
% 0   & 0      & \ddots & 0         \\
%     & 0      & 0      & B_N
% \end{array} \right)
% \] $$
%
% $head Hdia$$
% The argument $icode Hdia$$ is an $latex n \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex H_k \in \B{R}^{n \times n}$$ by
% $latex \[
%       H_k = Hdia(:, :, k)
% \] $$
%
%
% $head Hlow$$
% The argument $icode Hlow$$ is an $latex n \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex L_k \in \B{R}^{n \times n}$$ by
% $latex \[
%       L_k = Hlow(:, :, k)
% \] $$
%
% $head H$$
% The matrix $latex H$$ is defined by
% $latex \[
% H
% =
% \left( \begin{array}{cccc}
% H_1 & L_2^\R{T} & 0      &           \\
% L_2 & H_2    & L_3^\R{T} & 0         \\
% 0   & \ddots & \ddots & \ddots    \\
%     & 0      & L_N    & H_N
% \end{array} \right)
% \] $$
%
% $head ds$$
% The result $icode ds$$ is a column vector of length $latex r$$
% equal to the $latex \Delta s$$ components of the Newton step.
%
% $head dy$$
% The result $icode dy$$ is a column vector of length $latex p$$
% equal to the $latex \Delta y$$ components of the Newton step.
%
% $head du$$
% The result $icode du$$ is a column vector of length $latex r$$
% equal to the $latex \Delta u$$ components of the Newton step.
%
% $children#
%       example/newton_step_ok.m
% #$$
%
% $head Example$$
% The file $cref newton_step_ok.m$$ contains an example and test of
% $code ckbs_newton_step$$.
% It returns true if $code ckbs_newton_step$$ passes the test
% and false otherwise.
%
% $head Method$$
% The derivative of $latex F$$ is given by
% $latex \[
% F^{(1)} (s, y, u) =
% \left(
% \begin{array}{ccccc}
% D( 1_r ) &  0  & B  \\
% D( u )   & 0   & D(s) \\
% 0        & B^\R{T} & H
% \end{array}
% \right)
% \] $$
% It follows that
% $latex \[
% \left(\begin{array}{ccc}
% I & 0 & B \\
% U & S & 0\\
% 0 & B^T & C \\
% \end{array}\right)
% \left(\begin{array}{ccc}
% \Delta s \\ \Delta y \\ \Delta u
% \end{array}\right)
% =
% -
% \left(\begin{array}{ccc}
% s + b + By \\
% SU\B{1} - \mu\B{1}\\
% Cy + B^Tu + d
% \end{array}\right)\;.
% \] $$
% Below, $latex r_i$$ refers to row $latex i$$.
% Applying the row operations
% $latex \[
% \begin{array}{ccc}
% r_2 &=& r_2 - U r_1 \\
% r_3 &=& r_3 - B^T S^{-1} r_2
% \end{array}\;,
% \] $$
% we obtain the equivalent system
% $latex \[
% \left(\begin{array}{ccc}
% I & 0 & B \\
% 0 & S & -UB\\
% 0 & 0 & C + B^T S^{-1} U B
% \end{array}\right)
% \left(\begin{array}{ccc}
% \Delta s \\ \Delta y \\ \Delta u
% \end{array}\right)
% =
% -
% \left(\begin{array}{ccc}
% s + b + B y \\
% -U(b + B y) - \mu\B{1}\\
% Cy + B^T u + d + B^T S^{-1}
% \left(U(b + B y) + \mu \B{1}
% \right)
% \end{array}\right)\;.
% \] $$
% $latex \Delta y$$ is obtained
% from a single block tridiagonal solve (see third row of system).
% Then we immediately have
% $latex \[
% \Delta u = US^{-1}(b + B(y + \Delta y)) + \frac{\mu}{s}
% \] $$
% and
% $latex \[
% \Delta s = -s -b -B(y + \Delta y)\;.
% \] $$.
% $end
% ---------------------------------------------------------------------------
function [ds, dy, du] = ckbs_newton_step(mu, s, y, u, b, d, Bdia, Hdia, Hlow)
    m = size(Bdia, 1);
    n = size(Bdia, 2);
    N = size(Bdia, 3);
    p = n * N;
    r = m * N;
    %
    % check sizes
    if size(s,1) ~= r | size(u,1) ~= r | size(b,1) ~= r
        r
        size(s,1)
        size(u,1)
        size(b,1)
        error('ckbs_newton_step: argument sizes do not agree');
    end
    if size(y,1) ~= p | size(d,1) ~= p
        p
        size(y,1)
        size(d,1)
        error('ckbs_newton_step: argument sizes do not agree');
    end
    if size(Hdia,1) ~= n | size(Hdia,2) ~= n | ...
            size(Hlow,1) ~= n | size(Hlow,2) ~= n
        n
        size(Hdia,1)
        size(Hdia,2)
        size(Hlow,1)
        size(Hlow,2)
        error('ckbs_newton_step: argument sizes do not agree');
    end
    if size(Hdia,3) ~= N | size(Hlow,3) ~= N
        N
        size(Hdia,3)
        size(Hlow,3)
        error('ckbs_newton_step: argument sizes do not agree');
    end
    if min(s) <= 0
        error('ckbs_newton_step: min(s) <= 0 ');
    end
    if min(u) <= 0
        error('ckbs_newton_step: min(u) <= 0 ');
    end
    %
    % compute diagonal of D(u/s)
    u_sinv     = u ./ s;
    %
    % compute the diagonal blocks for H + B^T * D(u/s) * B
    Hplus = zeros(n, n, N);
    blk_m = 1 : m;
    for k = 1 : N
        Bk           = Bdia(:,:,k);
        Dk           = diag( u_sinv(blk_m) );
        Hplus(:,:,k) = Hdia(:,:,k) + Bk' * Dk * Bk;
        blk_m = blk_m + m;
    end
    %
    % compute right hand side in equation for Delta y
    Bt_u_d        = ckbs_blkdiag_mul_t(Bdia, u) + d;
    Cy            = ckbs_blktridiag_mul(Hdia, Hlow, y);
    By            = ckbs_blkdiag_mul(Bdia, y);
    row3          = Cy + Bt_u_d + ckbs_blkdiag_mul_t(Bdia, (u.*(b + By) + mu)./s);

    % compute du
    [dy, lambda]  = ckbs_tridiag_solve(Hplus, Hlow, -row3);

    Bydy = ckbs_blkdiag_mul(Bdia, y + dy);
    % compute du
    du            = (u.*(b + Bydy) + mu)./s;

    % compute ds
    ds            =  -s - b - Bydy;
    return
end
