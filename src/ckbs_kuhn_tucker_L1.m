% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_kuhn_tucker_L1$$ $newlinech %$$
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
%       L1
%       Aravkin
%       end end
% $$
%
% $section Compute Residual in Kuhn-Tucker Conditions for Robust L1$$
%
% $index ckbs_kuhn_tucker_L1$$
% $index kuhn_tucker_L1$$
%
% $index Kuhn-Tucker for L1, residual$$
% $index Tucker, Kuhn residual for L1$$
% $index residual, Kuhn-Tucker for L1$$
%
% $head Syntax$$
% $codei/[/F/] = ckbs_kuhn_tucker_L1(/
%       mu/, /s/, /y/, /r/, /b/, /d/, /Bdia/, /Hdia/, /Hlow/,
%       /pPlus/, /pMinus/)/$$
%
% $head Purpose$$
% This routine computes the residual in the Kuhn-Tucker conditions for
% the $latex \mu$$-relaxed affine L1 robust smoother problem.
%
% $head Problem$$
% Given
% $latex \mu \in \B{R}_+$$,
% $latex s \in \B{R}^{m \times N}$$,
% $latex y \in \B{R}^{n \times N}$$,
% $latex r \in \B{R}^{m \times N}$$,
% $latex b \in \B{R}^{m \times N}$$,
% $latex d \in \B{R}^{n \times N}$$,
% $latex B \in \B{R}^{m \times n \times N}$$,
% $latex H \in \B{R}^{n \times n \times N}$$,
% $latex p^+ \in \B{R}^{m \times N}$$,
% $latex p^- \in \B{R}^{m \times N}$$,
% the $latex \mu$$-relaxed affine L1 robust Kalman-Bucy smoother problem is:
%$latex \[
%\begin{array}{ll}
%{\rm minimize}
%& \frac{1}{2} y^\R{T} H(x) y + d(x)^\R{T} y
%       + \sqrt{\B{2}}^\R{T} (p^+ + p^-) -\mu \sum_{i =
%       1}^{mN} \log(p_i^+) - \sum_{i=1}^{mN} \mu \log(p_i^-)
%    \;{\rm w.r.t} \;  y \in \B{R}^{nN}\; , \; p^+ \in \B{R}_+^{M} \; , \; p^- \in \B{R}_+^{M}
%\\
%{\rm subject \; to} &  b(x) + B(x) y - p^+ + p^- = 0
%\end{array}
%\] $$
% In addition, $latex H$$ is symmetric block tri-diagonal with each block of
% size $latex n \times n$$ and
% $latex B$$ is block diagonal with each block of size $latex m \times n$$
%
% $head Lagrangian$$
% We use $latex r, \; s \in \B{R}^{m \times N}$$
% to denote $latex \mu /p^+\;,\; \mu/p^-$$, respectively, and
% we denote by $latex q$$ the lagrange multiplier associated to the
% equality constraint. We also use  $latex \B{1}$$
% to denote the vector of length $latex mN$$ with all its components
% equal to one, and $latex \B{\sqrt{2}}$$ to denote the vector of
% length $latex mN$$ with all its components equal to $latex
% \sqrt{2}$$.
% The corresponding Lagrangian is
% $latex \[
% L(p^+, p^-, y, q)  =
% \frac{1}{2} y^\R{T} H y + d^\R{T} y + \B{\sqrt{2}}^T(p^+ + p^-)
% - \mu \sum_{i=1}^{mN} \log(p_i^+) - \mu\sum_{i=1}^{mN}\log(p_i^-)
% + q^\R{T} (b + B y - p^+ + p^-)
% \] $$
% The partial gradients of the Lagrangian are given by
% $latex \[
% \begin{array}{rcl}
% \nabla_p^+ L(p^+, p^-, y, q ) & = & \B{\sqrt{2}} - q - r \\
% \nabla_p^- L(p^+, p^-, y, q) & = & \B{\sqrt{2}} + q - s \\
% \nabla_y L(p^+, p^-, y, q ) & = & H y + c + B^\R{T} q \\
% \nabla_q L(p^+, p^-, y, q ) & = & b + B y - p^+ + p^- \\
% \end{array}
% \] $$
% From the first two of the above equations,
% we have $latex q = (r - s)/2$$.
%
% $head Kuhn-Tucker Residual$$
% We use $latex D(s)$$
% to denote the diagonal matrix with $latex s$$
% along its diagonal.
% The Kuhn-Tucker Residual function
% $latex F : \B{R}^{4mN + nN} \rightarrow \B{R}^{4mN + nN}$$
% is defined by
% $latex \[
% F(p^+, p^-, r, s, y)
% =
% \left(
% \begin{array}{c}
% p^+ - p^- - b - B y       \\
% D(p^-) D(s) \B{1} - \tau \B{1} \\
% r + s - 2 \B{\sqrt{2}} \\
% D(p^+) D(r ) \B{1} - \tau \B{1}   \\
% H y + d + B^\R{T} (r - s)/2
% \end{array}
% \right)
% \] $$
% The Kuhn-Tucker conditions for a solution of the
% $latex \mu$$-relaxed constrained affine Kalman-Bucy smoother problem is
% $latex F(p^+, p^-, r, s, y) = 0 $$; see Equation (13) in
% $cref%Aravkin et al 2009%bib%Aravkin2009%$$
%
% $head mu$$
% The argument $icode mu$$ is a positive scalar specifying the
% relaxation parameter $latex \mu$$.
%
% $head s$$
% The argument $icode s$$ is an array of size $latex m \times
% N$$.
% All the elements of $icode s$$ are greater than zero.
%
% $head y$$
% The argument $icode y$$ is an array of size $latex n \times
% N$$
%
% $head r$$
% The argument $icode r$$ is an array of size $latex m \times
% N$$. All the elements of $icode r$$ are greater than zero.
%
% $head b$$
% The argument $icode b$$ is an array of size $latex m \times N$$.
%
% $head d$$
% The argument $icode d$$ is an array of size $latex n \times
% N$$.
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
% H_1 & L_2^\R{T} & 0         &           \\
% L_2 & H_2       & L_3^\R{T} & 0         \\
% 0   & \ddots    & \ddots    & \ddots    \\
%     & 0         & L_N       & H_N
% \end{array} \right)
% \] $$
%
% $head F$$
% The result $icode F$$ is a column vector of length $latex 4mN + nN$$
% containing the value of the
% $cref ckbs_kuhn_tucker_L1$$; Kuhn-Tucker L1 residual, i.e.,
% $latex F(p^+, p^-, s^+, s^-, y)$$.
%
% $children#
%       example/kuhn_tucker_L1_ok.m
% #$$
%
% $head Example$$
% The file $cref kuhn_tucker_L1_ok.m$$ contains an example and test of
% $code ckbs_kuhn_tucker_L1$$.
% It returns true if $code ckbs_kuhn_tucker_L1$$ passes the test
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------

function [F] = ckbs_kuhn_tucker_L1(mu, s, y, r, b, d, Bdia, Hdia, Hlow, pPlus, pMinus)

    % Assume the following about the input:

    % mu        scalar
    % s         m x N
    % y         n x N
    % r         m x N
    % b         m x N
    % d         n x N
    % Bdia     m x n x N
    % Hdia     n x n x N
    % Hlow  n x n x N
    % pPlus     m x N
    % pMinus    m x N

    %%
    % Determine the size of the problem
    n = size(y, 1);
    m = size(b, 1);
    N = size(y, 2);
    %

    %
    % check sizes
    if N ~= size(pPlus,2) || N ~= size(pMinus,2) || N ~= size(y,2) || ...
            N ~= size(r,2) || N ~= size(s,2) || N ~= size(d, 2) || ...
            N ~= size(Hdia, 3)|| N ~= size(Hlow, 3) || N ~= size(b, 2) ||...
            N ~= size(Bdia, 3)

        N

        size(pPlus,2)
        size(pMinus,2)
        size(y,2)
        size(r,2)
        size(s,2)
        size(d,2)
        size(Hdia,3)
        size(Hlow, 3)
        size(b, 2)
        size(Bdia, 3)


        error('L1_kuhn_tucker_L1: argument sizes do not agree');

    end
    if n ~= size(y,1) || n~= size(d, 1) || n ~= size(Hdia,1) || ...
            n ~= size(Hdia,2) || n ~= size(Hlow,1) || n ~= size(Hlow,2) || ...
            n ~= size(Bdia,2)
        n

        size(y,1)
        size(d,1)
        size(Hdia,1)
        size(Hdia,2)
        size(Hlow,1)
        size(Hlow,2)
        size(Bdia,2)

        error('L1_kuhn_tucker_L1: argument sizes do not agree');

    end
    if m ~= size(pPlus,1) || m ~= size(pMinus,1) || m ~= size(r,1) || m ~= size(s,1)||...
            m ~= size(b, 1) || m ~= size(Bdia, 1)
        m

        size(pPlus,1)
        size(pMinus,1)
        size(r,1)
        size(s,1)
        size(b, 1)
        size(Bdia, 1)

        error('L1_kuhn_tucker_L1: argument sizes do not agree');

    end




    yVec = y(:);
    dVec = d(:);
    pPlusVec   = pPlus(:);
    pMinusVec  = pMinus(:);
    rVec = r(:);
    sVec = s(:);

    bVec = b(:);

    % compute the necessary matrix products; i.e.,
    % H * y, B * y, and u^T * B
    Hy = ckbs_blktridiag_mul(Hdia, Hlow, y(:));
    Bt_SmR = ckbs_blkdiag_mul_t(Bdia, s(:)-r(:));
    By = ckbs_blkdiag_mul(Bdia, y(:));

    F = [ ...
        pPlusVec - pMinusVec - bVec - By; ...
        pMinusVec.*sVec - mu;...
        rVec + sVec - 2*sqrt(2); ...
        pPlusVec.*rVec - mu; ...
        Hy + dVec + 0.5*Bt_SmR; ...
        ];

    % F will be a vector of size (4m + n)*N
    if (4*m + n)*N ~= size(F, 1)
        (4*m + n)*N
        size(F, 1)
        error('L1_kuhn_tucker_L1: argument sizes do not agree');
    end
    if 1 ~= size(F, 2)
        size(F, 2)
        error('L1_kuhn_tucker_L1: argument sizes do not agree');
    end
end
