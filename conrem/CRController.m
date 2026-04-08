function CRMPC = CRController(ctrl,option)
% -------------------------------------------------------------------------
%                 [CRMPC] = CRController(ctrl,option)
% -------------------------------------------------------------------------
%  Function computes a Constraints Removal MPC policy based on given MPC
%  policy ctrl and selected options. For more details we refere users to 
%  read: https://github.com/holaza/mptplus/wiki
% -------------------------------------------------------------------------
% INPUTS:
%       ctrl              - MPT3 object (MPC controller)
% OPTIONS:
%       initializeType    -   CR controller tracks the last reduced MPC 
%                             problem. To recover all inequality cons, one 
%                             needs to call method ".InitializeCRMPC". 
%                             This parameter changes the was how the CR 
%                             controller performs this inicialization:
%                               1: Initialize with only feasible ineq. cons
%                               0: Initialize with all ineq. cons (also 
%                                  infeasible)
%       CRfixed           -   Parameter allows user to enable/disable 
%                             constraints removal. (This is helpful when 
%                             measuring evaluation time of the MPC method 
%                             with/withoud CR removal.)
%                               1: CR MPC does not remove any inequality cons
%                               0: CR MPC does remove any inequality cons
%       SolveType         -   Option how to time the (quadprog)
%                             optimisation problem.
%                               1: timeit()
%                               0: tic-toc (or for the higher Matlab 
%                               versions the build-in quadprog solver time)       
%       SigmaShiftBias    -   Shifts the computed sigma values by: sigma + Bias
%       SigmaShiftGain    -   Shifts the computed sigma values by: sigma * gain
%                             (The formula: sigma_new = sigma*Gain + Bias)
% OUTPUTS:
%      CRMPC                 - the Constrains Removal MPC policy (MPT3 object)
% -------------------------------------------------------------------------
% % EXAMPLE:
% model = LTISystem(ss([1 1; 0 1], [1; 0.5], [1 0], 0,'Ts', 0.1));
% model.x.max = [5; 5];
% model.x.min = [-5; -5];
% model.u.max = 1;
% model.u.min = -1;
% model.x.penalty = QuadFunction(diag([1 1]));
% model.u.penalty = QuadFunction(0.1);
% model.x.with('terminalPenalty');
% model.x.terminalPenalty = model.LQRPenalty();
% model.x.with('terminalSet');
% model.x.terminalSet = model.LQRSet();
% 
% ctrl = MPCController(model,5);
% crmpc = CRController(ctrl,{'SigmaShiftBias',10});
% 
% data1 = ctrl.simulate(x0, Nsim);
% data2 = crmpc.simulate(x0, Nsim);
% 
% figure
% subplot(3,1,1)
% hold on
% plot(data1.X')
% plot(data2.X', '--r')
% subplot(3,1,2)
% hold on
% stairs(data1.U)
% stairs(data2.U, '--r')
% subplot(3,1,3)
% stairs(data2.ncons)
% -------------------------------------------------------------------------

if nargin < 1, error('MPTplus: Ups, not enough inputs were selected'); end

validScalarNonegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
isBinary = @(x) ismember(x,[0 1]);

ip = inputParser;
ip.addParameter('initializeType', 1, isBinary);
ip.addParameter('CRfixed', 0, isBinary);
ip.addParameter('SolveType', 1, isBinary);
ip.addParameter('SigmaShiftBias', 0, validScalarNonegNum);
ip.addParameter('SigmaShiftGain', 1, validScalarNonegNum);
if nargin == 2, ip.parse(option{:}); else ip.parse(); end
options = ip.Results;

%% Check compatibility of MPC formulations 
% The user provides an MPT3 MPC policy with a specific MPC formulations 
% (e.g. delta-u, y-tracking, ...). This part of hte code validates of the
% given formulations are compatible with the CR approach.

MPCFormulation = DetectMPCFormulations(ctrl);

formulation = [];
if MPCFormulation.deltau
    formulation = 'delta-u';
elseif MPCFormulation.xtracking
    formulation = 'xtracking';
elseif MPCFormulation.ytracking
    formulation = 'ytracking';
end
if ~isempty(formulation)
    error(['MPTplus: We are sorry, but the ', ...
            formulation, ...
            ' formulation is not ' ...
            'supported with the Constraints removal techniques.'])
end

%% Recover optimization problem from the ctrl object
% Reformulate the ctrl problem as:
% min  z'*Q*z + b'*z
% s.t. Aineq * z <= bineq
%        Aeq * z == beq
%      Ainit * z == z0 
% where z = [{y_1,..y_N}, yref; {u_1,...,u_N}, uprev, x_0, {x_1,...,x_N}, xref]
% and z0 = [x0, uprev, xref, yref] is vector of parameters
fprintf('Processing the optimisation problem ...\n')
optprob = RedefineOptimizationProblem(ctrl);

%% Declare Lyapunov function
% TODO: For future releases create also a helper Lyapunov function 
%      (that will be evaluated right after the state measurement)

% Case 1: Lyapunov = obj function
L = @(z) z'*optprob.Q*z + optprob.b'*z;

%% Assign sigma values
if options.CRfixed
    fprintf('Assigning sigma values is not needed (CRfixed = 1)\n')
    sigma = inf(optprob.nineq,1);
    sigmaCompInfo = [];
else
    fprintf('Assigning sigma values based on the Lyapunov function ...\n')
    [sigma, sigmaCompInfo] = AssignSigmaValues(optprob,L);
    fprintf('\t ... done in %.2f seconds. \n',sigmaCompInfo.execution_time)

    % Shift sigmavalues
    sigma = sigma/options.SigmaShiftGain;
    sigma = sigma - options.SigmaShiftBias;
end

%% Sorting constraints based on their sigma values (ascending order)
if ~options.CRfixed
    [sigma,indx] = sort(sigma,'ascend');
    optprob.Aineq = optprob.Aineq(indx,:);
    optprob.bineq = optprob.bineq(indx,:);
    optprob.ineq = @(z) optprob.Aineq*z <= optprob.bineq;
    optprob.const = @(z,x0)[optprob.eq(z);optprob.ineq(z);optprob.initcon(z,x0)];
end

%% Create CR MPC policy
secOutputs.options = options;
secOutputs.optprob = optprob;
secOutputs.L = L;
secOutputs.sigma = sigma;
secOutputs.sigmaCompInfo = sigmaCompInfo;

CRMPC = CRMPCController(ctrl,secOutputs);

end











function [optprob] = RedefineOptimizationProblem(ctrl)
% ctrl -> MPC policy (MPT3 object)

% if ctrl.wasModified() || isempty(ctrl.optimizer)
%     % Re-generate the optimizer if the problem setup was changed
%     ctrl.construct();
% end

% Just to create the optimizer
ctrl.construct();
internal_model = ctrl.optimizer.model;
% exctract the original optimization problem:
% min  z'*Q*z + b'*z
% s.t. A(i,2:end) * z <= A(i,1)     % if A(i, 1) ~= 0
%      A(i,2:end) * z == 0          % if A(i, 1) == 0
% where z = [x,u]
A = internal_model.F_struc;
b = internal_model.c;
Q = internal_model.Q;

% Reformulate the problem as:
% min  z'*Q*z + b'*z
% s.t. Aineq * z <= bineq
%        Aeq * z == beq
%      Ainit * z == z0 
% where z = [{y_1,..y_N}, yref; {u_1,...,u_N}, uprev, x_0, {x_1,...,x_N}, xref]
% and z0 = [x0, uprev, xref, yref] is vector of parameters
optProbDef = ['min  z''*Q*z + b''*z; s.t.: Aineq * z <= bineq; ', ...
    'Aeq * z == beq; Ainit * z == z0; ', ...
    'where z = [{y_1,..y_N}, yref; {u_1,...,u_N}, uprev, x_0, {x_1,...,x_N}, xref]',...
    ' and z0 = [x0, uprev, xref, yref] is vector of parameters. ' ...
    '(particular order -> see optprob.requestedFormat.indices)'];
indx_eq = find(A(:,1) == 0);
indx_in = find(A(:,1) ~= 0);
Aeq = A(indx_eq,2:end);
beq = A(indx_eq,1);
Aineq = A(indx_in,2:end);
bineq = A(indx_in,1);
% Define parameters as Ainit*z == x0
thetaIndices = internal_model.parameterIndex;
Ainit = zeros(length(thetaIndices),size(Aeq,2));
for k = 1:length(thetaIndices)
    Ainit(k,thetaIndices(k)) = 1;
end

% Define useful function 
symbvar = @() sdpvar(size(Q,1),1);
cost = @(z) z'*Q*z + b'*z;
ineq = @(z) Aineq*z <= bineq;
eq = @(z) Aeq*z == beq;
initcon = @(z,theta) z(thetaIndices) == theta;
% initcon = @(z,x0) Ainit*z == x0;
const = @(z,theta) [eq(z); ineq(z); initcon(z,theta)];



% Further analysis of the optimisation problem
MPCFormulation = DetectMPCFormulations(ctrl);
optinfo.symbvar = symbvar;
optinfo.cost = cost;
optinfo.const = const;
requestedFormat = Prepare_MPT_requestedFormat_xyu(ctrl,MPCFormulation,optinfo);
variables = Define_variables_xyu(ctrl,requestedFormat);

% Saving results
optprob.optProbDef = optProbDef;
optprob.Q = Q;
optprob.b = b;
optprob.Aineq = Aineq;
optprob.bineq = bineq;
optprob.Aeq = Aeq;
optprob.beq = beq;
optprob.Ainit = Ainit;
optprob.symbvar = symbvar;  % sdpvar(size(Q,1),1)
optprob.cost = cost;        % z'*Q*z + b'*z
optprob.ineq = ineq;        % Aineq*z <= bineq
optprob.eq = eq;            % Aeq*z == beq
optprob.initcon = initcon;  % z(internal_model.parameterIndex) == x0
optprob.const = const;      % [eq(z); ineq(z); initcon(z,x0)]
optprob.neq = size(optprob.Aeq,1);
optprob.nineq = size(optprob.Aineq,1);
optprob.nc = optprob.neq + optprob.nineq;
optprob.nx = ctrl.nx;
optprob.nu = ctrl.nu;
if isprop(ctrl.model,'C')
    optprob.ny = ctrl.model.ny;
else
    optprob.ny = [];
end
optprob.nz = size(optprob.Aeq,2);
optprob.requestedFormat = requestedFormat;
optprob.MPCFormulation = MPCFormulation;
optprob.variables = variables;
end



function [requestedFormat] = Prepare_MPT_requestedFormat_xyu(ctrl,MPCFormulation,optinfo)

% Creating requestedFormat for future evaluation
try 
    requestedFormat = ctrl.Internal.yalmipData.internal.requestedFormat;
catch
    if isprop(ctrl.model,'C')
        names = {'cost' 'u'  'x'  'y'};
        components = {'model'  'u'  'x'  'y'};
        dims = {[1 1]  [ctrl.model.nu ctrl.N]  [ctrl.model.nx ctrl.N+1]  [ctrl.model.ny ctrl.N]};
    else
        names = {'cost' 'u'  'x'};
        components = {'model'  'u'  'x'};
        dims = {[1 1]  [ctrl.model.nu ctrl.N]  [ctrl.model.nx ctrl.N+1]};
    end
    requestedFormat.names = names;
    requestedFormat.names_upper = ([names{1}, names(2:end)]);
    requestedFormat.components = components;
    requestedFormat.dims = dims;
end


% ------------------------------- indices----------------------------------
fprintf('Identifying indices of the QP ...\n')
[indices,theta] = getIndices(ctrl, MPCFormulation, optinfo);
fprintf('\t... done!\n')
%   x0 - indices
model = ctrl.optimizer.model;
requestedFormat.indices = indices;
% requestedFormat.indices.x0 = model.parameterIndex + 1;
requestedFormat.indices.theta = theta;
% ------------------------------------------------------------------------


% ----------------------------- ZINDICES ---------------------------------
requestedFormat.zindices.Vector = requestedFormat.names(2:end);
requestedFormat.zindices.Dims = requestedFormat.dims(2:end);
requestedFormat.zindices.x = requestedFormat.indices.x - 1;
if isfield(requestedFormat.indices,'y')
    requestedFormat.zindices.y = requestedFormat.indices.y - 1;
else
    requestedFormat.zindices.y = [];
end
requestedFormat.zindices.u = requestedFormat.indices.u - 1;

% requestedFormat.zindices.x0 = model.parameterIndex(1:ctrl.nx)';
if MPCFormulation.xtracking
    requestedFormat.zindices.xref = requestedFormat.indices.theta.xref - 1;
    requestedFormat.zindices.Vector{length(requestedFormat.zindices.Vector)+1} = 'xref';
    requestedFormat.zindices.Dims{length(requestedFormat.zindices.Dims)+1} = [ctrl.nx 1];
else
    requestedFormat.zindices.xref = nan;
end
if MPCFormulation.ytracking
    requestedFormat.zindices.yref = requestedFormat.indices.theta.yref - 1;
    requestedFormat.zindices.Vector{length(requestedFormat.zindices.Vector)+1} = 'yref';
    requestedFormat.zindices.Dims{length(requestedFormat.zindices.Dims)+1} = [ctrl.ny 1];
else
    requestedFormat.zindices.yref = nan;
end
if MPCFormulation.deltau
    requestedFormat.zindices.uprev = requestedFormat.indices.theta.uprev;
    requestedFormat.zindices.Vector{length(requestedFormat.zindices.Vector)+1} = 'uprev';
    requestedFormat.zindices.Dims{length(requestedFormat.zindices.Dims)+1} = [ctrl.nu 1];
else
    requestedFormat.zindices.uprev = nan;
end
requestedFormat.zindices.theta = model.parameterIndex;
% ------------------------------------------------------------------------

% ------------------------------ Theta -----------------------------------
% Create also indices w.r.t. vector of parameters
requestedFormat.thetaindices.x0 = [];
requestedFormat.thetaindices.xref = [];
requestedFormat.thetaindices.yref = [];
requestedFormat.thetaindices.uprev = [];

% X0
Vector = {'x0'};
Dims = {[ctrl.nx 1]};
requestedFormat.thetaindices.x0 = requestedFormat.indices.theta.x0 - 1;
% Delta-u 
if MPCFormulation.deltau
    requestedFormat.thetaindices.uprev = requestedFormat.indices.theta.uprev - 1;
    Vector{length(Vector)+1} = 'uprev';
    Dims{length(Dims)+1} = [ctrl.nu 1];
end
% X-tracking
if MPCFormulation.xtracking
    requestedFormat.thetaindices.xref = requestedFormat.indices.theta.xref - 1;
    Vector{length(Vector)+1} = 'xref';
    Dims{length(Dims)+1} = [ctrl.nx 1];
end
% Y-tracking
if MPCFormulation.ytracking
    requestedFormat.thetaindices.yref = requestedFormat.indices.theta.yref - 1;
    Vector{length(Vector)+1} = 'yref';
    Dims{length(Dims)+1} = [ctrl.ny 1];
end
requestedFormat.thetaindices.Vector = Vector;
requestedFormat.thetaindices.Dims = Dims;

internal_model = ctrl.optimizer.model;
thetaIndices = internal_model.parameterIndex;
requestedFormat.thetaindices.theta = thetaIndices;
requestedFormat.thetaindices.thetainfo = 'theta = [x0, uprev, xref, yref]';
end



function [variables] = Define_variables_xyu(ctrl,requestedFormat)
variables.name = requestedFormat.zindices.Vector;
variables.dim.nu = length(requestedFormat.zindices.u);
variables.dim.ny = length(requestedFormat.zindices.y);
variables.dim.nx = length(requestedFormat.zindices.x);
variables.dim.theta = length(requestedFormat.zindices.theta) ;
temp = variables.dim.nu + variables.dim.ny + variables.dim.nx + variables.dim.theta;
variables.dim.nz = temp - ctrl.nx;
end



function [sigma,sigmaCompInfo] = AssignSigmaValues(optprob,L)
% AssignSigmaValues computes sigma values for each inequality constraint
%
% sigma(k) = minimal value of the objective if k-th inequality is active
%
% Inputs:
%   optprob.Aineq   - QP inequality matrix (Aineq*x <= bineq)
%   optprob.bineq   - QP inequality vector
%   optprob.Aeq     - QP equality matrix (Aeq*x == beq)
%   optprob.beq     - QP equality vector
%   optprob.variables.dim.nz - number of decision variables
%   L(z)            - function handle for objective
%
% Output:
%   sigma           - vector of sigma values (length = number of inequalities)

% number of constraints
nc = size(optprob.Aineq,1);

% Preallocate sigma
sigma = nan(nc,1);

% initialize results of each QP solved
resinfo.infeas = 0;
resinfo.suboptimal = 0;
resinfo.numericalIssues = 0;

% Loop over all inequalities
tStart = tic;
for k = 1:nc
    [sigma(k),resinfo] = computeSigmaForKQuadprog(optprob,k,resinfo);
    ShowProgress('Calculations:     ', k, nc)
end
execution_time = toc(tStart);
sigmaCompInfo.execution_time = execution_time;
resinfo.feasible =  nc - resinfo.infeas - resinfo.suboptimal - resinfo.numericalIssues;
sigmaCompInfo.resinfo = resinfo;
end



function [sigma,resinfo] = computeSigmaForKQuadprog(optprob,k,resinfo)
% ================================================================
% QUADRATIC PROGRAM solved by MATLAB quadprog
%
% The solver computes the vector z that minimizes
%
%     1/2 * z' * H * z + f' * z
%
% subject to the constraints
%
%     A * z <= b                (linear inequality constraints)
%     Aeq * z == beq            (linear equality constraints)
%     lb <= z <= ub             (lower and upper bounds)
%
% where
%
%     z     ... decision variable (n x 1 vector)
%     H     ... Hessian matrix (n x n, symmetric)
%     f     ... linear cost term (n x 1)
%     A     ... inequality constraint matrix
%     b     ... inequality constraint vector
%     Aeq   ... equality constraint matrix
%     beq   ... equality constraint vector
%     lb    ... lower bounds (optional)
%     ub    ... upper bounds (optional)
%
% In this implementation the problem
%
%     min  z' * Q * z + b' * z
%
% is converted to quadprog form using
%
%     H = 2*Q
%     f = b
%
% ================================================================

% Activate k-th inequality as equality
Aeq = [optprob.Aeq; optprob.Aineq(k,:)];
beq = [optprob.beq; optprob.bineq(k)];

% Remove k-th inequality from Aineq
mask = true(size(optprob.Aineq,1),1);
mask(k) = false;
Aineq = optprob.Aineq(mask,:);
bineq = optprob.bineq(mask);

% ===== COST FUNCTION =====
H = 2*optprob.Q;          % quadprog uses 1/2*z'*H*z + f'*z
f = optprob.b;

% ===== CONSTRAINTS =====

% inequality constraints
A = Aineq;
b_qp = bineq;

% equality constraints
Aeq_qp = Aeq;
beq_qp = beq;

% ===== INITIAL GUESS (needed for active-set) =====
z_init = zeros(size(optprob.Q,1),1);   % replace with a feasible guess if required

% ===== OPTIONS =====
if exist('mosekopt.m','file')
    options = [];
else
    options = optimoptions('quadprog');
    options.Algorithm = 'active-set';
    options.Display = 'off';
end

% ===== SOLVE =====
[~, fval, exitflag, ~] = quadprog( ...
    H, f, ...          % cost
    A, b_qp, ...       % inequality
    Aeq_qp, beq_qp, ...% equality
    [], [], ...        % bounds (none here)
    z_init, ...        % initial guess
    options);

%  1	Optimal solution found 
%  0	Maximum number of iterations reached (solution not guaranteed optimal)
% -2	Problem is infeasible (constraints cannot be satisfied)
% -3	Problem is unbounded (objective can decrease indefinitely)
% Robust handling
if exitflag ~= 1
    if exitflag == -2 || exitflag == -3
        sigma = inf;  
        resinfo.infeas = resinfo.infeas + 1;
    elseif exitflag == 0
        sigma = -inf; % suboptimal
        fprintf('Constraint %d: suboptimal solution -> this con may not be removed.\n',k);
        resinfo.suboptimal = resinfo.suboptimal + 1;
    else
        error('MPTplus: ups quadporg returns unidentified exitflag.')
    end
else
    sigma = fval;
end

end



function [MPCFormulation] = DetectMPCFormulations(ctrl)
% Functionidentifies identifies formulations of MPC givne by the user.
%
% Inputs:
%       ctrl            - MPT3 MPC controller
%
% Outputs:
%   MPCFormulation.deltau       - delta-u formulation
%   MPCFormulation.xtracking    - x-tracking formulation
%   MPCFormulation.ytracking    - y-tracking formulation
%   MPCFormulation.xinitFormat  - format of inintial conditions
% -------------------------------------------------------------------------

Y = ctrl.toYALMIP.internal.xinitFormat;
MPCFormulation.deltau = any(strcmp(Y.names,'u.previous'));
MPCFormulation.xtracking = any(strcmp(Y.names,'x.reference'));
MPCFormulation.ytracking = any(strcmp(Y.names,'y.reference'));
MPCFormulation.xinitFormat = Y;
end



function [indices,theta] = getIndices(ctrl, MPCFormulation, optinfo)
% Functionidentifies indices of the optimized vector Z (of the exported QP)
% with respect to the MPC formulation of MPT3.
%
% Inputs:
%       ctrl            - MPT3 MPC controller
%       MPCFormulation  - mpc formulations (given in the ctrl)
%       optinfo.symbvar - vector of optimized varialbes (yalmip sdpvars)
%       optinfo.cost    - objective function of the QP (yalmip formulation)
%       optinfo.const   - constraints of the QP (yalmip formulation)
%
% Outputs:
%       indices    - openloop sequence (containing X/U/Y)
%       theta       - vector of inputs to the MPT3 MPC object
% -------------------------------------------------------------------------

for m = 1:5
    [openloop,Z,~,~] = GenerateASolution(ctrl,MPCFormulation,optinfo);
    [indices] = findIndicesInYopt_tol(openloop, Z);
    % Check if we have found indices (if not try it again)
    hasEmpty = 0;
    f = fieldnames(indices);
    for k = 1:length(f)
        if isempty(indices.(f{k}))
            hasEmpty = 1;
            break;
        end
    end
    % If indices are correctly identified then continue
    if hasEmpty, break; end
end
[theta] = identifyParameters(ctrl, indices);
end



function [openloop,Z,theta,thetaz] = GenerateASolution(ctrl,MPCFormulation,optinfo)
% Function generates initial condition for then given MPC 'ctrl'. Then both
% the MPC is solved returning the openloop sequence. Finally the associated
% QP problem is solved to return the optimized vector Z.
%
% Inputs:
%       ctrl            - MPT3 MPC controller
%       MPCFormulation  - mpc formulations (given in the ctrl)
%       optinfo.symbvar - vector of optimized varialbes (yalmip sdpvars)
%       optinfo.cost    - objective function of the QP (yalmip formulation)
%       optinfo.const   - constraints of the QP (yalmip formulation)
%
% Outputs:
%       openloop    - openloop sequence (containing X/U/Y)
%       Z           - vector of otimized variables (containing X/U/Y/theta)
%       theta       - vector of inputs to the MPT3 MPC object
%       thetaz      - vector of input parameters to the QP
% -------------------------------------------------------------------------

x0 = [];
yref = [];
xref = [];
uprev = [];
for k = 1:50
    x0 = zeros(ctrl.nx,1) + (rand(ctrl.nx,1).^2)*sign(rand-0.5)/2;
    theta = {x0};
    if MPCFormulation.ytracking
        yref = zeros(ctrl.ny,1) + (rand(ctrl.ny,1).^2)*sign(rand-0.5)/2;
        yref = min(max(ctrl.model.y.min, yref), ctrl.model.y.max);
        theta = [theta, {'y.reference', yref}];
    end
    if MPCFormulation.xtracking
        xref = zeros(ctrl.nx,1) + (rand(ctrl.nx,1).^2)*sign(rand-0.5)/2;
        xref = min(max(ctrl.model.x.min, xref), xref);
        theta = [theta, {'x.reference', xref}];
    end
    if MPCFormulation.deltau
        uprev = zeros(ctrl.nu,1) + (rand(ctrl.nu,1).^2)*sign(rand-0.5)/2;
        uprev = min(max(ctrl.model.u.min, uprev), ctrl.model.u.max);
        theta = [theta, {'u.previous', uprev}];
    end

    % Solve the MPC for the given initial condition
    [~,feas,openloop] = ctrl.evaluate(theta{:});

    if feas
        % TODO: double check this !!!
        % that a for the optimisation problem is: [x0; uprev; xref; yref]
        thetaz = [x0; uprev; xref; yref];

        % Solve the QP problem for the given initial condition
        opt_settings = sdpsettings('verbose',0);
        z = optinfo.symbvar();
        obj = optinfo.cost(z);
        con = optinfo.const(z,thetaz);
        yopt = optimize(con,obj,opt_settings);
        Z = double(z);
        if ~yopt.problem
            % both prolbems are feasible
             break;
        end
    end
    if k == 50 
        error('MPTplus: can not choose a random feasible condition, try to run me again ...')
    end
end
end



function [indices] = findIndicesInYopt_tol(openloop, Z)
% This function identifies vectors of X/Y/U from the optimized vector Z.
% The approach is to simply compare results from the openloop with Z.
%
% Inputs:
%       openloop.U   - vector of inputs
%       openloop.X   - vector of states
%       openloop.Y   - vector of outputs
%       Z            - vector of otimized variables (containing X/U/Y/theta)
%
% Outputs:
%       indices.J   - indices of the objective function
%       indices.u   - indices of inputs
%       indices.x   - indices of states
%       indices.y   - indices of outputs
% -------------------------------------------------------------------------

% Identifying indexing of z (vector of optimized variables)
tol = 1e-8;

f = lower(fieldnames(openloop));
% Remove fields not stored in Yopt
f(strcmp(f,'runtime')) = [];
f(strcmp(f,'cost')) = [];

idx = struct('x', [], 'y', [], 'u', []);
for k = 1:length(f)

    name = f{k};
    value = openloop.(upper(name));
    v = value(:);           % convert to vector
    n = length(v);

    found = false;

    for i = 1:(length(Z)-n+1)

        block = Z(i:i+n-1);

        if all(abs(block - v) < tol)
            idx.(name) = i:i+n-1;
            found = true;
            break
        end
    end

    if ~found
        idx.(name) = [];
        warning(['No match found for field ', name])
    end
end

% Shift it by +1 as we are considering indices of the MPT3 object
indices.J = 1;
f = fields(idx);
for k = 1:length(f)
    indices.(f{k}) = idx.(f{k}) + 1;
end
end



function [thetas] = identifyParameters(ctrl, idx)
% This function identifies vectors of parameters theta from the optimized 
% vector Z. Simple heuristic is used: yref is near vector of Y, uprev is
% near vector of U, ...
%
% Inputs:
%       ctrl    - MPT3 MPC controller
%       idx.u   - indices of inputs
%       idx.x   - indices of states
%       idx.y   - indices of outputs
% 
% Outputs:
%       thetas.yref   - indices of y-references
%       thetas.uprev  - indices of previous contorl actions
%       thetas.x0     - indices of initial states
%       thetas.xref   - indices ofx-references
% -------------------------------------------------------------------------

% predefine the parameters
theta = struct( ...
    'yref', [], ...
    'uprev', [], ...
    'x0', [], ...
    'xref', [] );

% What do we know is indices of input parameters:
model = ctrl.optimizer.model;
p = model.parameterIndex(:);

% ---- number of states from requestedFormat ----
nx = ctrl.model.nx;

% ---- beginning of X block ----
firstX = min(idx.x-1);

% ---- x0 = first nx indices at beginning of X ----
theta.x0 = p(p >= firstX & p < firstX + nx);

% remove x0 from parameter list
p(ismember(p,theta.x0)) = [];

% ---- remaining parameters ----
for i = 1:length(p)

    ind = p(i);

    % next to Y -> yref
    if any(abs(ind - idx.y-1) == 1)
        theta.yref(end+1) = ind;

        % next to U -> uprev
    elseif any(abs(ind - idx.u-1) == 1)
        theta.uprev(end+1) = ind;

        % next to X but not x0 -> xref
    elseif any(abs(ind - idx.x-1) == 1)
        theta.xref(end+1) = ind;
    end
end

% Shift it by +1 as we are considering indices of the MPT3 object
thetas = struct( ...
    'yref', [], ...
    'uprev', [], ...
    'x0', [], ...
    'xref', [] );
f = fields(theta);
for k = 1:length(f)
    thetas.(f{k}) = theta.(f{k}) + 1;
end
end



function ShowProgress(operation, k, kf)
% Inputs:
%          operation - sting
%          k         - current iteration
%          kf        - final number of iterations
%          ttime     - time of the previous operation

% if k == 1, fprintf(1, strcat(operation,'       :')); end
 if k == 1, fprintf(1, strcat(operation,'    :')); end

kkf = fix(k/kf*100);
if kkf < 10
    fprintf(1, ' \b\b\b\b%d %%', kkf);
else
    fprintf(1, ' \b\b\b\b\b%d %%', kkf);
end

if k == kf
    fprintf('\n');
end
end



