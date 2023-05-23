classdef PolynomialMPC < EMPCController
    % Class representing explicit polynomial-approximation-based controller
    %
    % Constructor:
    %   ctrl = PolynomialMPC(explicit_controller,options)
    %
    %   options
    %       .aprx_degree - approximation degree of the resulting polynomial
    %                      (Default value: 3)
    %       .polya_degree - degree of the polya polynomial (Default value: 1)
    %       .gamma        - enforced Lyapunov decay rate (default gamma = 1e-5)
    %                       for regions around the origin:
    %                       V(x(k+1)) - V(x(k)) <= gamma * ||x(k)||_p
    %       .norm_p       - norm p (= 1 or inf) used in the enforced Lyapunov
    %                       decay rate (default p = inf)
    %       .dV_eps       - enforced Lyapunov decay rate (default dV_eps = 1e-6)
    %                       for regions not around the origin:
    %                       V(x(k+1)) - V(x(k)) <= -dV_eps
    %       .merge        - =0 (default) no intermediate merging of the XUset
    %                       =1 simple intermediate greedy merging
    %                       =2 intermediate optimal merging
    %       .convex       - =0 don't try to find a convex representation
    %                          of polytope arrays in intermediate steps
    %                       =1 (default) find a convex representation of
    %                          for polytope arrays in intermediate steps
    %       .psolvers     - which projection method to use
    %                       (default is [3 1 5 2 4 0], see help
    %                       polytope/projection for more details)
    %       .rel_tol      - relative tollerance used e.g. to verify a set
    %                       difference (Default 1e-6)
    %
    % Literature: M. Kvasnica – J. Löfberg – M. Fikar: Stabilizing polynomial
    % approximation of explicit MPC. Automatica, vol. 10 (47), pp. 2292–2297,
    % 2011.

    % Current limitations:
    %   LTI system is assumed: i.e. no ULTI or hybrid dynamics
    %   MPT3 creates PWQ obj function oif terminal set is set --- WHAT ???


    properties
        polynomial
        controller
        XUset
        info
    end

    methods

        function out = getName(obj)
            out = 'Explicit polynomial-approximation-based controller';
        end

        function obj = PolynomialMPC(ctrl, option)
            % Constructor:
            %   clip = PolynomialMPC(explicit_controller, option)

            if nargin == 0
                return
            elseif ~isa(ctrl, 'EMPCController')
                error('Input must be an EMPCController object.');
            end

            % set options/default values
            ip = inputParser;
            ip.addParameter('aprx_degree', 3);
            ip.addParameter('polya_degree', 1);
            ip.addParameter('gamma', 1e-5);
            ip.addParameter('dV_eps', 1e-6);
            ip.addParameter('merge', 0);
            ip.addParameter('convex', 1);
            ip.addParameter('norm_p', inf);
            ip.addParameter('verbose', 1);
            ip.addParameter('rel_tol', 1e-6);
            ip.addParameter('psolvers', [3 1 5 2 4 0]);
            if nargin == 2, ip.parse(option{:}); else ip.parse(); end
            options = ip.Results;

            % copy data from the source object
            obj.importUserData(ctrl.copy());

            % Construct polynomial-approximation-based MPC Controller
            % compute the XU stabilizing set
            XUset = mpt_XUctrl_mpt3(ctrl, options);
            % degree of the approximation polynomial
            aproximation_degree = options.aprx_degree;
            % degree of the polya polynomial
            polya_degree = options.polya_degree;
            % compute the approximation
            [polynomial, controller, info] = xu_approx_robust_mpt3(ctrl, XUset, aproximation_degree, polya_degree);

            % Assign PROPERTIES
            obj.polynomial = polynomial;
            obj.controller = controller;
            obj.XUset = XUset;
            obj.info = info;
        end

        function [u, feasible, openloop] = evaluate(obj, xinit, varargin)
            % Evaluates the clipping controller for a given point
            %
            % u = controller.evaluate(x0) evaluates the explicit MPC
            % solution and returns the approximated (suboptimal) control input
            % associated to point "x0". If "x0" is outside of the
            % controller's domain, "u" will be NaN.
            %
            % [u, feasible] = controller.evaluate(x0) also returns a
            % feasibility boolean flag.
            %
            %% TODO: Check if the following works correctly:
            %
            % [u, feasible, openloop] = controller.evaluate(x0) also
            % returns the full open-loop optimizer in "openloop.U" and the
            % optimal cost in "openloop.cost". Moreover, "openloop.X" and
            % "openloop.Y" will be set to NaN. This is due to the fact that
            % storing open-loop predictions of states and outputs in the
            % explicit solution would considerably increase its size.
            %
            % u = controller.evaluate(x0, 'x.reference', ref, 'u.prev', u0)
            % includes "ref" and "u0" into the vector of initial
            % conditions.

            % evaluate the approximated control law
            u = xu_poly_eval(obj.polynomial, xinit);

            % Feasibility check
            if( obj.partition.Domain.contains(xinit) == 0 )
                u = NaN;
                feasible = 0;
            else
                feasible = 1;
            end

            openloop.cost = NaN;
            openloop.U = NaN;
            openloop.X = NaN;
            openloop.Y = NaN;
            openloop.info = 'OPENLOOP evaluation of polynomial-based approximated explicit MPC is not supported, yet.';

        end

        function plot(obj,plotSaturated)
            % plot the approximated solution
            if ( nargin < 1 )
                error('Not enought inputs')
            end
            if ( nargin < 2 )
                plotSaturated = 0;
            end
            xu_plot_mpt3(obj,obj.polynomial,obj.XUset,plotSaturated)
        end
    end

end



function [p_opt, ctrl, runtime] = xu_approx_robust_mpt3(ctrl, XUset, order, polya_degree, gridreg)
% find polynomial u(x)=a'*x s.t. [x; u] \in XU(i) \forall x \in P(i), \forall i

if nargin < 4, polya_degree = 1; end
if nargin < 5, gridreg = 10; end
usecost = true;
yalmip('clear');

XU = XUset.XU;

nx = ctrl.nx;
a = sdpvar(ones(1, order+1), repmat(nx, 1, order+1));
a{end} = 0;

E = cell(1,ctrl.nr);
for i = 1:ctrl.nr
    E{i} = ctrl.partition.Set(i).V;
end

objrob = 0;
Frob = set([]);

fprintf('robustifying constraints...\n');
tic
for i = 1:ctrl.nr
    obj = 0;
    F = [];
    w = sdpvar(size(E{i}, 1), 1);
    F = F + [uncertain(w)];
    F = F + [sum(w)==1; w >= 0];
    x = 0;
    for j = 1:size(E{i}, 1)
        x = x + E{i}(j, :)'*w(j);
    end
    u = xu_poly_eval(a, x);
    F = F + [ismember([x; u], XU(i))];
    x0 = grid(ctrl.partition.Set(i), gridreg);
    for ii = 1:size(x0, 1)
        u_opt = ctrl.evaluate(x0(ii, :)');
        u_aprx = xu_poly_eval(a, x0(ii, :)');
        obj = obj + norm(u_opt - u_aprx, 1);
    end
    [FF, OO] = robustify(F, obj, sdpsettings('robust.polya', polya_degree));
    Frob = Frob + FF;
    objrob = objrob + OO;
end
if isa(objrob, 'double') || ~usecost
    objrob = [];
end
runtime.robustify = toc;
fprintf('\nsolving the problem...\n');
tic
sol = optimize(Frob, objrob, sdpsettings('robust.polya', polya_degree));
runtime.lp = toc;
runtime.sol = sol;

% extract optimal values of the parameters of the approximation polynomial
p_opt = cell(1, length(a));
for i = 1:length(a)
    p_opt{i} = double(a{i});
end
end


function [out] = mpt_XUctrl_mpt3(ctrl,Options,LyapFcn)
% MPT_XUCTRL computes the stabilizing control set on the basis of a given
% PWA Lyapunov function
%
% XUset = mpt_XUsetPWALyap(ctrl,Options);
% XUset = mpt_XUsetPWALyap(ctrl,Options,LyapFcn); [NOT IMPLEMENTED YET!!]
%
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% computes the stabilizing controller set U(x) on the basis of a given PWA
% Lyapunov function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl                    - Explicit controller (MPTCTRL object)
%
% LyapFcn                 - PWA Lyapunov Function (if not provided: the value
%                           function in ctrl is assumed to be Lyapunov for the
%                           controlled system)
%                           NOTE: The LyapFcn needs to be defined over the
%                                 same polyhedral partition as ctrl.Pn!!
%                                 (e.g. MPT_LYAPUNOV creates such a LyapFcn)
%  LyapFcn.Pn             - Polyhedral partition of the state-space
%  LyapFcn.Bi, LyapFcn.Ci - cell arrays containing the PWA Lyapunov function
%
% Options.gamma           - enforced Lyapunov decay rate (default gamma = 1e-5)
%                           for regions around the origin:
%                            V(x(k+1)) - V(x(k)) <= gamma * ||x(k)||_p
% Options.norm_p          - norm p (= 1 or inf) used in the enforced Lyapunov
%                           decay rate (default p = inf)
% Options.dV_eps          - enforced Lyapunov decay rate (default dV_eps = 1e-6)
%                           for regions not around the origin:
%                            V(x(k+1)) - V(x(k)) <= -dV_eps
% Options.merge           - =0 (default) no intermediate merging of the XUset
%                           =1 simple intermediate greedy merging
%                           =2 intermediate optimal merging
% Options.convex          - =0 don't try to find a convex representation
%                              of polytope arrays in intermediate steps
%                           =1 (default) find a convex representation of
%                              for polytope arrays in intermediate steps
% Options.psolvers        - which projection method to use
%                           (default is [3 1 5 2 4 0], see help
%                           polytope/projection for more details)
%
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% out                     - XU sets details:
%   .details.XU           - XU sets
%   .details.XUproj       - XU sets projected down to X-space
%   .details.XUtime       - time spent when computing XU sets
%   .details.proj_time    - time spent when computing projection of XU sets
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
%
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the
%          Free Software Foundation, Inc.,
%          59 Temple Place, Suite 330,
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

% Set Default values
if nargin < 2, Options = []; end
if ~isfield(Options, 'psolvers'), Options.psolvers = [3 1 5 2 4 0]; end
if nargin < 3, LyapFcn = []; end

% Create the XU set
startt = clock;
XU = mpt_XUsetPWALyap_mpt3(ctrl, Options, LyapFcn);
xutime = etime(clock, startt);

% Check boundness of the XU set:
if ~isBounded(XU)
    error([mfilename ': XU set is not bounded, probably due to numerical problems.']);
end
startt = clock;
XUproj = projection(XU, 1:ctrl.model.nx, Options);
if ~isBounded(XUproj)
    error([mfilename ': projected XU set is not bounded, probably due to numerical problems.']);
end
projectiontime = etime(clock, startt);

% Define outputs
out.XU = XU;
out.XUproj = XUproj;
out.XUtime = xutime;
out.proj_time = projectiontime;
end




function XUset = mpt_XUsetPWALyap_mpt3(ctrl,Options,LyapFcn)
% MPT_XUSETPWALYAP computes the stabilizing control set on the basis of
% a given PWA Lyapunov function.
%
% XUset = mpt_XUsetPWALyap(ctrl,Options);
% XUset = mpt_XUsetPWALyap(ctrl,Options,LyapFcn); [NOT IMPLEMENTED YET!!]
%
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% computes the stabilizing controller set U(x) on the basis of a given PWA
% Lyapunov function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl                    - Explicit controller (MPTCTRL object)
%
% LyapFcn                 - PWA Lyapunov Function, i.e. PolyUnion object with
%                           a PWA function 'obj' defined over it.
%                           See 'ctrl.cost' polyunion as an example (if not
%                           provided: the value function in ctrl is assumed
%                           to be Lyapunov for the controlled system)
%                           NOTE: The LyapFcn needs to be defined over the
%                                 same polyhedral partition !!(e.g.
%                                 MPT_LYAPUNOV from mpt2 creates such a LyapFcn)
%
% Options.gamma           - enforced Lyapunov decay rate (default gamma = 1e-5)
%                           for regions around the origin:
%                            V(x(k+1)) - V(x(k)) <= gamma * ||x(k)||_p
% Options.norm_p          - norm p (= 1 or inf) used in the enforced Lyapunov
%                           decay rate (default p = inf)
% Options.dV_eps          - enforced Lyapunov decay rate (default dV_eps = 1e-6)
%                           for regions not around the origin:
%                            V(x(k+1)) - V(x(k)) <= -dV_eps
% Options.merge           - =0 (default) no intermediate merging of the XUset
%                           =1 simple intermediate greedy merging
%                           =2 intermediate optimal merging
% Options.convex          - =0 don't try to find a convex representation
%                              of polytope arrays in intermediate steps
%                           =1 (default) find a convex representation of
%                              for polytope arrays in intermediate steps
% Options.psolvers        - which projection method to use
%                           (default is [3 1 5 2 4 0], see help
%                           polytope/projection for more details)
% Options.rel_tol         - relative tolerance (default 1e-6)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% XUset                   - Polyhedral partition of the state+input-space,
%                           where [x' U(x)']' \in XUset

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
%
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the
%          Free Software Foundation, Inc.,
%          59 Temple Place, Suite 330,
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

% check intputs to this function
narginchk(1,3);
if ~strcmp(ctrl.getName,'Explicit MPC controller')
    error('Only explicit controllers are supported by this function.');
end

% ----------------------------- !!! TODO !!! -----------------------------
% Choosing Lyap function (Currently, we use only the value function.)
if ~isempty(LyapFcn)
    warning(['Custom LyapFcn was provided! Make sure that the ' ...
        'lyapunov function is defined over',...
        ' the same polyhedral partion as the controller partition',...
        ' (with the same ordering of the polytopes)! Moreover, ' ...
        'the function name has to be "obj"!'])
else
    LyapFcn = ctrl.cost;
end

% Forcing ndyn = 1 (i.e. only LTI systems are assumed / no hybrid dynamics)
% forcing nbool = 0 (Try to determine from the mpt3 contorller if it contains any bool variables)
ndyn = 1;
nbool = 0;

if nbool > 0
    error('Systems with boolean inputs are not supported.');
end

if ndyn > 1
    error('Sorry, PWA systems are not supported')
end

% if ctrl.model.x.penalty.type == 2 % ZZZ
if ( isfield(ctrl.model.x.penalty, 'H') )
    error('Polynomial approximation does not support 2-norm solutions, yet.');
end
% -------------------------------------------------------------------------


% detect flat regions, because
%  if Pn(ii) is lowdim => XUset is lowdim  (not considered in the soution)
%  if Pn(jj) is lowdim => XUset is lowdim  (not considered in the soution)
% !!! TODO: include if U-set is also flat => same reasoning!!!
ind_flat = [];
for ii = 1:ctrl.nr
    xout = ctrl.partition.Set(ii).chebyCenter;
    Rn = xout.r;
    if Rn < Options.rel_tol
        Bset = ctrl.partition.Set(ii).outerApprox();
        d = abs(Bset.Internal.ub-Bset.Internal.lb);
        if any(d > Options.rel_tol)
            ind_flat = [ind_flat ii];
        end
    end
end
if ~isempty(ind_flat)
    fprintf('%d flat regions will be skipped.\n', length(ind_flat));
end

% Determine indices of regions that contains the origin
% (due to different handling during the main loop)
[~, zero_inside, ~] = contains(ctrl.partition,zeros(ctrl.nx,1));
zero_inside = sort(zero_inside);

% extreme values of the LyapFunction, i.e. currently J(x), over all regions
Jmaxmin = nan(ctrl.nr,2);
for ii=1:ctrl.nr
    H = LyapFcn.Set(ii).A;
    K = LyapFcn.Set(ii).b;
    Bi = LyapFcn.Set(ii).getFunction('obj').F;
    Ci = LyapFcn.Set(ii).getFunction('obj').g;
    [~,Jmin] = mpt_solveLPs(Bi,H,K,[],[],[],[]);
    [~,Jmax] = mpt_solveLPs(-Bi,H,K,[],[],[],[]);
    Jmaxmin(ii,:) = [-Jmax Jmin] + Ci;
end

isConvexOpt.abs_tol = Options.rel_tol;
if isfield(Options, 'usehull')
    isConvexOpt.usehull = Options.usehull;
end

% Initialize the final set of polyhedra as []
XUset = [];

% main loops
LoopInputs.ctrl = ctrl;
LoopInputs.LyapFcn = LyapFcn;
LoopInputs.ind_flat = ind_flat;
LoopInputs.Options = Options;
LoopInputs.isConvexOpt = isConvexOpt;
if ndyn < 2 % LTI system
    for ii = 1:ctrl.nr  % source region
        if ismember(ii,ind_flat), continue; end % skip "flat" regions

        XUsetIter = [];
        jj_ind = find(Jmaxmin(ii,1) > Jmaxmin(:,2));

        if ismember(ii,zero_inside)
            % Case 1: (0 \in Pi) & (Jj_min < Ji_max)
            %         => XUset might be non-convex
            IterCase = 1;
            XUsetIter_case1 = MainLoopCases(LoopInputs,IterCase,jj_ind,ii,Jmaxmin(ii,2));
            % Store results of this itteration
            XUsetIter = [XUsetIter XUsetIter_case1];

        else % CASE 2 & 3: ~(0 \in Pi)
            % CASE 2:  ~(0 \in Pi) & (0<= Jj(x) <= Ji_min)
            %       => XUset is convex!!
            jj_case2 = jj_ind(Jmaxmin(jj_ind,2)<Jmaxmin(ii,2));
            IterCase = 2;
            XUsetIter_case2 = MainLoopCases(LoopInputs,IterCase,jj_case2,ii,Jmaxmin(ii,2));

            % CASE 3:  ~(0 \in Pi) & (Ji_min < Jj(x) <= Ji_max)
            %       => XUset is possibly non-convex!!
            jj_case3 = setdiff(jj_ind,jj_ind(Jmaxmin(jj_ind,1)<=Jmaxmin(ii,2)));
            IterCase = 3;
            XUsetIter_case3 = MainLoopCases(LoopInputs,IterCase,jj_case3,ii,Jmaxmin(ii,2));

            % Post processing of CASE2 & CASE3
            XUsetIter = [XUsetIter_case2 XUsetIter_case3];
            % if set is convex, use it
            if Options.convex
                U = PolyUnion(XUsetIter);
                if isConvex_mod(U,isConvexOpt.abs_tol)
                    XUsetIter = U.convexHull;
                end
            end
        end % different cases for ii

        % Intermediate Merging: post-merging for region ii
        if Options.merge>0 & length(XUsetIter)>1
            if Options.merge ==1      % greedy merging
                XUsetIter = merge(XUsetIter);
            elseif Options.merge ==2  % optimal mergin
                XUsetIter = merge(XUsetIter,struct('greedy',0));
            end
        end
        XUset = [XUset XUsetIter];
    end
else  % PWA system
    error('Sorry, PWA systems are not supported ...')
end
end






% SUB-FUNCTIONS __________________________________
function P = XUpoly(Hu, Ku, H1, K1, H2, K2, b1, c1, b2, c2, A, B, f,...
    nx, nu, Jmin, Jmax, dV_eps, Gx, Gu, Gc)

H = [];
K = [];

nc1 = length(K1);
nc2 = length(K2);

% x \in P1
H = [H1 zeros(nc1, nu)];
K = K1;

% A*x + B*u + f \in P2
H = [H; H2*A H2*B];
K = [K; K2 - H2*f];

% u \in U
H = [H; zeros(2*nu, nx) Hu];
K = [K; Ku];

if nargin > 18
    % PWA guards:
    % Gx*x + Gu*u <= Gc
    H = [H; Gx Gu];
    K = [K; Gc];
end

if Jmin > 0   % for Case 3: 0 < Ji_min <= J(x^+)
    H = [H; -b2*[A B]];
    K = [K; c2 - Jmin - b2*f];
end

if isinf(Jmax)
    % b2*(A*x+B*u+f) + c2 - (b1*x + c1) <= - dV_eps
    H = [H; (b2*A - b1) b2*B];
    K = [K; c1 - c2 - dV_eps - b2*f];

else
    % b2*(A*x+B*u+f) + c2 - Jmax <= - dV_eps
    H = [H; b2*A b2*B];
    K = [K; -c2 + Jmax - dV_eps - b2*f];
end

P = Polyhedron(H, K);
end


% (0 \in Pi): regions around the origin
function P = XUpoly_zero(Hu, Ku, H1, K1, H2, K2, b1, c1, b2, c2, A, B, f,...
    nx, nu, gamma, p, psolvers, Gx, Gu, Gc)

H = [];
K = [];

nc1 = length(K1);
nc2 = length(K2);

% x \in P1
H = [H1 zeros(nc1, nu)];
K = K1;

% A*x + B*u + f \in P2
H = [H; H2*A H2*B];
K = [K; K2 - H2*f];

if nargin > 18
    % PWA guards:
    % Gx*x + Gu*u <= Gc
    H = [H; Gx Gu];
    K = [K; Gc];
end

% u \in U
H = [H; zeros(2*nu, nx) Hu];
K = [K; Ku];

% b2*(A*x+B*u+f) + c2 - (b1*x + c1) <= - gamma * ||x||_p
if isinf(p)  % inf-norm case
    H = [H zeros(size(H,1),1)];

    H = [H; (b2*A - b1) b2*B gamma];
    K = [K; c1 - c2 - b2*f];

    H = [H; eye(nx) zeros(nx,nu) -ones(nx,1); -eye(nx) zeros(nx,nu) -ones(nx,1)];
    K = [K; zeros(2*nx,1)];

else % 1-norm case
    H = [H zeros(size(H,1),nx)];

    H = [H; (b2*A - b1) b2*B gamma*ones(1,nx)];
    K = [K; c1 - c2 - b2*f];

    H = [H; eye(nx) zeros(nx,nu) -eye(nx); -eye(nx) zeros(nx,nu) -eye(nx)];
    K = [K; zeros(2*nx,1)];
end

P = Polyhedron(H, K);
if isFullDim(P)
    for ii = 1:length(psolvers)
        Pproj = projection(P,[1:(nx+nu)], [], struct('psolvers', psolvers(ii:end)));
        if ~isFullDim(Pproj) | ~isBounded(Pproj)
            if length(psolvers) > ii
                fprintf('Warning: Numerical problems with projection, recomputing...\n');
            else
                error('Numerical problems with projection.');
            end
        else
            % hopefully correct (or at least fully-dimensional) projection found
            break
        end
    end
    P = Pproj;
end
end



function ts = isConvex_mod(obj,convxtol)
%
% check if the polyhedron array forms a convex union
%

if nargin < 2, convxtol = 1e-6;end

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isConvex;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
end

% if isempty(obj.Internal.Convex)
% compute the convex hull
H = obj.convexHull();
% the union is convex if H\U = emptyset
ts = all(isEmptySet(mldivide(H, obj.Set, true)));

% This part is new: Need to better impose tollerance
if ~ts
    temp = mldivide(H, obj.Set, true);
    temp2 = temp.chebyCenter;
    if temp2.r <= convxtol, ts = 1;end
end
return

end % function



function xu_plot_mpt3(ctrl,aopt,XUset,plotSaturated)
% plot the approximated solution
if nargin < 3, error('Not enought inputs'); end
if nargin < 4, plotSaturated = 0; end

order = length(aopt)-1;
umax = ctrl.model.u.max;
umin = ctrl.model.u.min;

if ctrl.nx == 1
    x = ctrl.partition.Domain.grid(300);
else
    error('Currently only 1D partitions can be plotted.');
end
y_sat = zeros(1, size(x, 1));
y = zeros(1, size(x, 1));
z = zeros(1, size(x, 1));
for j = 1:size(x, 1)
    y_sat(j) = xu_poly_eval(aopt, x(j, :)', umin, umax);
    y(j) = xu_poly_eval(aopt, x(j, :)');
    % evaluate the original PWA control law
    u_open = ctrl.feedback.getFunction('primal').feval(x(j, :)');
    z(j) = u_open(1:ctrl.model.nu);
end

figure
hold on
XUset.XU.plot('color','y','linewidth', 2)
plot(x, z, 'b--', 'LineWidth', 3);
if plotSaturated
    plot(x, y, 'g', 'LineWidth', 3);
    plot(x, y_sat, 'r--', 'LineWidth', 3);
else
    plot(x, y, 'r--', 'LineWidth', 3);
end
axis([min(x)*1.1 max(x)*1.1 min(y_sat)*1.1 max(y_sat)*1.1])
if plotSaturated, s = '(saturated)'; else, s = ''; end
title(sprintf('Approximation order: %d %s', order, s));
hold off
end


function p = xu_poly_eval(a, x, umin, umax)
order = length(a) - 1;
p = sum(a{end});
for i = 1:order
    p = p + a{i}*x.^(order+1-i);
end
if nargin > 3
    %     p = max(min(p, umax), umin);
    p = 0.5*(0.5*(p + umax - abs(p - umax)) + umin + abs(0.5*(p + umax - abs(p - umax)) - umin));
end
end


function [XUsetIter_case] = MainLoopCases(LoopInputs,IterCase,jj_ind,indx,Jmaxmin)
% Extract common inputs
ctrl = LoopInputs.ctrl;
LyapFcn = LoopInputs.LyapFcn;
ind_flat = LoopInputs.ind_flat;
Options = LoopInputs.Options;
isConvexOpt = LoopInputs.isConvexOpt;

% Lyap values:
H1 = LyapFcn.Set(indx).A;
K1 = LyapFcn.Set(indx).b;
b1 = LyapFcn.Set(indx).getFunction('obj').F;
c1 = LyapFcn.Set(indx).getFunction('obj').g;

% !!! We are not assuming PWA dynamics !!!
A = ctrl.model.A;
B = ctrl.model.B;
f = ctrl.model.f;

% System dimensions
nx = ctrl.nx;
nu = ctrl.nu;

% constraints on U
if isprop(ctrl.model.u,'setConstraint')
    Hu = ctrl.model.u.setConstraint.A;
    Ku = ctrl.model.u.setConstraint.b;
else
    Hu = [eye(nu); -eye(nu)];
    Ku = [ctrl.model.u.max; -ctrl.model.u.min];
end

% verbouse
if Options.verbose > 0
    if IterCase == 1
        if Options.gamma > 0
            disp(sprintf('region: i=%d /%d  [region at 0: additional projection computations]',...
                indx,ctrl.nr))
        else
            disp(sprintf('region: i=%d /%d  [region at 0]',indx,ctrl.nr))
        end
    elseif IterCase == 2
        disp(sprintf('region: i=%d /%d',indx,ctrl.nr))
    end
end

XUsetIter_case = [];
for jj = jj_ind'
    if ismember(jj,ind_flat), continue; end % skip "flat" regions
    H2 = LyapFcn.Set(jj).A;
    K2 = LyapFcn.Set(jj).b;
    b2 = LyapFcn.Set(jj).getFunction('obj').F;
    c2 = LyapFcn.Set(jj).getFunction('obj').g;

    switch IterCase
        case 1
            if Options.gamma > 0
                P = XUpoly_zero(Hu, Ku, H1, K1, H2, K2, b1, c1, b2, c2, A, B, f,...
                    nx, nu, Options.gamma, Options.norm_p, Options.psolvers);
            else % don't use projection
                P = XUpoly(Hu, Ku, H1, K1, H2, K2, b1, c1, b2, c2, A, B, f,...
                    nx, nu, 0, inf, 0);
            end
        case 2
            P = XUpoly(Hu, Ku, H1, K1, H2, K2, b1, c1, b2, c2, A, B, f,...
                nx, nu, 0, Jmaxmin, Options.dV_eps);
        case 3
            P = XUpoly(Hu, Ku, H1, K1, H2, K2, b1, c1, b2, c2, A, B, f,...
                nx, nu, Jmaxmin, inf, Options.dV_eps);
    end

    XUsetIter_case = [XUsetIter_case P];
end

% if set is convex, use it (does not need to be convex for IterCase = {1, 3}!)
if Options.convex
    U = PolyUnion(XUsetIter_case);
    if isConvex_mod(U,isConvexOpt.abs_tol)
        XUsetIter_case = U.convexHull;
    elseif IterCase == 2
        error('Set should be convex (numerical problems)')
    end
end

end
