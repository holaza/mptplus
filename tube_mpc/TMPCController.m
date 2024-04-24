function TMPC = TMPCController(model,N,option)
% -------------------------------------------------------------------------
%                 [TMPC] = TMPCController(model,N,option)
% -------------------------------------------------------------------------
%  Function computes a Tube MPC policy based on given uncertain model,
%  prediction horizon, and selected options. For more details we refere
%  users to read:
%  https://github.com/holaza/mptplus/wiki
% ------------------------------------------------------------------------
% INPUTS:
%      model                - MPT3 object with ULTISystem
%      N                    - Prediction horizon
% OPTIONS:
%      Tube_tol           - allowed tollerace (default = 1e-4 )
%      Tube_MaxIter       - maximal number of iterations (default = 1e3)
%      solType              - {0/1} defines the retuned optimized values
%                             (the default value is 1):
%                             0: TMPC/TEMPC returns exact TMPC outputs:
%                                           [u_opt, x_opt]
%                             1: TMPC/TEMPC returns compact TMPC outputs:
%                                     Utube = u_opt + K*( x - x_opt ),
%                                where [u_opt, x_opt] are the original TMPC
%                                outputs, K is the LQR for the given
%                                model (given as u = K*x), and x is the
%                                state measurement.
%     LQRstability          - {0} Does not compute any terminal penalty/set
%                             {1} Computes the terminal penalty/set w.r.t.
%                                 the LTI system (not the ULTI system)
%     TubeType              - {'implicit', 'explicit'} Type of the tube 
%                             embedded in the optimization problem. By
%                             default, TubeType = 'explicit'.
%     Nz                    - number of terminal steps. By default, Nz = 0.
%                             (i.e. prediction steps with control actions 
%                              computed as u_k = Kz*x, where Kz is also 
%                              given by the user)
%     Kz                    - terminal LQR 
%     Pz                    - terminal cost (when Kz is used) 
%
% OUTPUTS:
%      TMPC                 - tube MPC policy (MPT3 object)
% ------------------------------------------------------------------------
% % EXAMPLE:
% % LTI system
% model = ULTISystem('A', [1, 1; 0, 1], 'B', [0.5; 1], 'E', [1, 0; 0, 1]);
% % Input constraints
% model.u.min = [-1];
% model.u.max = [ 1];
% % State constraints
% model.x.min = [-100; -100];
% model.x.max = [ 100;  100];
% % Disturbance constraints
% model.d.min = [-0.1; -0.1];
% model.d.max = [ 0.1;  0.1];
% % Penalty functions
% model.x.penalty = QuadFunction(diag([1, 1]));
% model.u.penalty = QuadFunction(diag([0.01]));
% % Prediction horizon
% N = 9;
% % Options: Tube in the implicit form + terminal penalty and constraints
% option = {'TubeType','implicit','LQRstability',1}; 
% % Compute the Tube MPC object
% TMPC = TMPCController(model,N,option);
% % Perform simulation
% x0 = [-5;-2]; % initial condition
% Nsim = 15;    % number of simulation steps
% data = TMPC.simulate(x0,Nsim);
% figure
% plot(data.X')
% % Computation of Explicit solution
% eMPC = TMPC.toExplicit();
% u_implicit = eMPC.evaluate(x0)  % Explicit MPC evaluation
% [u, feasible, openloop] = eMPC.evaluate(x0) % Enriched output
% ------------------------------------------------------------------------

%% Checks / initial setup
% Perform checks:
if nargin < 2, error('MPT+: Ups, not enough inputs were selected'); end
% We do not support parametric uncertainties
if iscell(model.A) || iscell(model.B)
    error('MPT+: Implemented Tube MPC design method does not support parametric uncertainties (LPV systems), yet.')
end

% Define option validation handles
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
isBinary = @(x) ismember(x,[0 1]);
isValidMatrix = @(x) all([all(~isnan(x(:))),all(isnumeric(x(:)))]);
isValidTubeType = @(x) any(strcmpi(x,{'explicit','implicit'}));

% set options/default values
ip = inputParser;
% ip.addParameter('TEMPC', 0);
ip.addParameter('adjList', 0, isBinary);
ip.addParameter('Tube_tol', 1e-4, validScalarPosNum);
ip.addParameter('Tube_MaxIter', 1e3, validScalarPosNum);
ip.addParameter('solType', 1, isBinary);
ip.addParameter('LQRstability', 0, isBinary);
ip.addParameter('TubeType', 'explicit', isValidTubeType);
ip.addParameter('Nz', 0, validScalarPosNum);
ip.addParameter('Kz', nan, isValidMatrix);
ip.addParameter('Pz', nan, isValidMatrix);
if nargin == 3, ip.parse(option{:}); else ip.parse(); end
options = ip.Results;

%% TMPC setup
% copy model to do not alter the original object
workingModel = model.copy(); 
% Construct common sets and LQR for both implicit/explicit tube approaches
secOutputs = ConstructSetsAndLQR(workingModel);
% Store the presiction horizon
secOutputs.N = N;
% Store and verify provided option inputs
secOutputs = PreprocessOptions(workingModel,secOutputs,options);
% Prepare all ingredients for the implicit (and also explicit) tube MPC
secOutputs = PreProcessing_Common(workingModel,secOutputs,options);
% Precompute all necessary components for implicit/explicit tube approaches
if strcmpi(options.TubeType,'implicit')
    secOutputs = PreProcessing_ImplicitTube(workingModel,secOutputs);
    secOutputs = DefineMPCProblem_ImplicitTube(secOutputs,workingModel,options);
elseif strcmpi(options.TubeType,'explicit')
    secOutputs = PreProcessing_ExplicitTube(workingModel,secOutputs);
    secOutputs = DefineMPCProblem_ExplicitTube(secOutputs,workingModel,options);
end

%% Construction of MPC
% Preconsturct an MPC object
if options.solType
    TMPC = MPCController(workingModel,1);
else
    TMPC = MPCController(workingModel,N);
end
TMPC.construct();

% Construct the Tube (implicit) MPC
TMPC = TubeMPCController( TMPC , secOutputs );

% Overwrite the optimizer based on the solType
if options.solType
    TMPC.optimizer = secOutputs.OptimProblem.MPCoptimizers.compact;
else
    TMPC.optimizer = secOutputs.OptimProblem.MPCoptimizers.exact;
end

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [secOutputs] = ConstructSetsAndLQR(model)
% Function creates sets (constraints) and LQR for further processing.
%
% INPUTS:
%           model                 - model of the system (as MPT3 object)
% OUTPUTS:
%           secOutputs.K          - nominal LQR policy (as u = K*x)
%           secOutputs.Pw         - Bounded additive uncertainty
%           secOutputs.Px         - State constraints
%           secOutputs.Pu         - Input constraints
%           secOutputs.Pxu        - State-Input constraints
%           secOutputs.Pdu        - Delta-u constratins

% Computing LQR for the nominal model (not the ULTI!) as u = K*x
K = -dlqr(model.A,model.B,model.x.penalty.weight,model.u.penalty.weight); 
secOutputs.K = K; 

% Bounded additive uncertainty
if isprop(model.d,'setConstraint')
    W = model.d.setConstraint;
else
    W = Polyhedron('lb',model.d.min,'ub',model.d.max);
end
W = model.E*W;
Wset = Polyhedron('lb',-max(abs(W.V)),'ub',max(abs(W.V)));
secOutputs.Pw = NormalizePolyhedron(Wset);

% State constraints
if isprop(model.x,'setConstraint')
    Px = model.x.setConstraint;
else
    Px = Polyhedron('lb',model.x.min,'ub',model.x.max);
end
secOutputs.Px = NormalizePolyhedron(Px);

% Input constraints
if isprop(model.u,'setConstraint')
    Pu = model.u.setConstraint;
else
    Pu = Polyhedron('lb',model.u.min,'ub',model.u.max);
end
secOutputs.Pu = NormalizePolyhedron(Pu);

% State-input constraints
secOutputs.Pxu = secOutputs.Px*secOutputs.Pu;

% Delta-u constraints
if all([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
    dUset_full = Polyhedron('lb',model.u.deltaMin ,'ub',model.u.deltaMax );
    secOutputs.Pdu = NormalizePolyhedron(dUset_full);
elseif isprop(model.u,'deltaMin')
    dUset_full = Polyhedron('lb',model.u.deltaMin );
    secOutputs.Pdu = NormalizePolyhedron(dUset_full);
elseif isprop(model.u,'deltaMax')
    dUset_full = Polyhedron('ub',model.u.deltaMax );
    secOutputs.Pdu = NormalizePolyhedron(dUset_full);
else
    secOutputs.Pdu = {};
end
end



function [PP] = NormalizePolyhedron(P)
if any(P.b <= 0), error('MPT+: Ups, normalization of the polyhedron failed!'); end
PP = Polyhedron('A',P.A.*P.b.^(-1),'b',P.b.*P.b.^(-1));
end



function [secOutputs] = PreprocessOptions(model,secOutputs,options)
secOutputs.Nz = options.Nz;
secOutputs.Kz = options.Kz;
secOutputs.Pz = options.Pz;

% Is Nz valid
if any([isnan(options.Nz),~isnumeric(options.Nz)])
    error('MPT+: Parameter "Nz" must be numeric')
elseif options.Nz > secOutputs.N
    error('MPT+: Parameter "Nz" cannot be larger than the prediciton horizon "N"!')
elseif options.Nz == secOutputs.N
    warning('MPT+: Parameter "Nz" is equal to the prediction horizon "N".')
end


% Check Kz, Pz inputs only if Nz > 0
if options.Nz > 0
    % Both Kz and Pz need to be given
    if any([sum(isnan(options.Kz(:))),sum(isnan(options.Pz(:)))])
        issuestr = 'MPT+: not enough inputs were selected. When requesting Nz > 0, then it is also needed to provide: \n';
        if isnan(options.Kz)
            issuestr = [issuestr, '\t\tLQR terminal gain "Kz" \n'];
        end
        if isnan(options.Pz)
            issuestr = [issuestr, '\t\tLQR terminal penalty "Pz" \n'];
        end
        error(issuestr)
    end

    % Are Kz/Pz of right dimensions?
    [t1,t2] = size(options.Kz);
    if any([t1 ~= model.nu,t2 ~= model.nx])
        error('MPT+: Provided matrix "Kz" has invalid dimensions!')
    end
    [t1,t2] = size(options.Pz);
    if any([t1 ~= t2,t2 ~= model.nx])
        error('MPT+: Provided matrix "Pz" has invalid dimensions!')
    end

    % Is the Kz stabilizing?
    if any(abs(eig(model.A + model.B*options.Kz)) >= 1)
        warning('MPT+: Non-stabilizing LQR gain "Kz" was provided, i.e., (abs(eig(A+B*Kz))<1) does not hold.')
    end

    % Is the Pz positive semi-definite?
    if any(eig(options.Pz) < 0)
        warning('MPT+: Provided terminal penalty "Pz" is not semi-positive define, i.e., all(eig(options.Pz) <= 0) does not hold.')
    end
end
end



function [secOutputs] = PreProcessing_Common(model,secOutputs,options)
% Function preprocess variables for further processing.
%
% INPUTS:
%           model                   - model of the system (as MPT3 object)
%           secOutputs.Pw           - Bounded additive uncertainty
%           secOutputs.Pxu          - State-Input constraints
%           options.Tube_tol        - tollerence for the (implicit) MRPIS
%           options.Tube_MaxIter    - maximal number of iterations for MRPIS
% OUTPUTS:
%           secOutputs.Ns           - number of iterations to determine MRPIS 
%           secOutputs.alpha        - computed tollerence for the MRPIS
%           secOutputs.alpha_tol    - final tolerance for ALPHA-convergence
%           secOutputs.f            - implicit tube

[Ns, alpha, alpha_tol ] = alpha_Ns(model,secOutputs,options.Tube_tol,options.Tube_MaxIter);

f = support_f(model,secOutputs,Ns,alpha);

secOutputs.Ns = Ns;
secOutputs.alpha = alpha;
secOutputs.alpha_tol = alpha_tol;
secOutputs.f = f;
secOutputs.nu = model.nu; % This is needed only for the evlauation function
secOutputs.nx = model.nx; % This is needed only for the evlauation function
end



function [secOutputs] = PreProcessing_ImplicitTube(model,secOutputs)
% Function preprocess variables for further processing.
%
% INPUTS:
%           model                   - model of the system (as MPT3 object)
%           secOutputs.Px           - Bounded additive uncertainty
%           secOutputs.Pu           - State-Input constraints
%           secOutputs.f            - implicit tube
% OUTPUTS:
%           secOutputs.Tube         - minimal robust positive invariant set
%           secOutputs.Px_robust    - robustustified state constraint
%           secOutputs.Pu_robust    - robustustified input constraint
%           secOutputs.Pdu_robust   - robustustified delta input constraint

% Reminder:
% Implicit tube 'secOutputs.f' represents how much we need to cut from 
% state-input sets [x,u] (in that particular order).

fnx = size(secOutputs.Px.H,1);
fnu = size(secOutputs.Pu.H,1);
fx = secOutputs.f(1:fnx);
fu = secOutputs.f(fnx+1:fnx+fnu);

% Robust positive invariant set
Tube = {};

% State constraints robustification
Px_robust = Polyhedron('A',secOutputs.Px.A,'b',secOutputs.Px.b - fx);
if Px_robust.isEmptySet
    error('MPT+: Cannot continue as the state constraint set is empty!')
end
% Input constraints robustification
Pu_robust = Polyhedron('A',secOutputs.Pu.A,'b',secOutputs.Pu.b - fu);
if Pu_robust.isEmptySet
    error('MPT+: Cannot continue as the input constraint set is empty!')
end
% Delta input constraints robustification
if any([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
    Pdu_robust = Polyhedron('A',secOutputs.Pdu.A,'b',secOutputs.Pdu.b - fu);
end

% Store all computed sets
secOutputs.Tube = Tube;
secOutputs.Px_robust = Px_robust;
secOutputs.Pu_robust = Pu_robust;
if any([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
    secOutputs.Pdu_robust = Pdu_robust;
end
end




function [secOutputs] = PreProcessing_ExplicitTube(model,secOutputs)
% Function preprocess variables for further processing.
%
% INPUTS:
%           model                   - model of the system (as MPT3 object)
%           secOutputs.Px           - State constraints
%           secOutputs.Pu           - Input constraints
%           secOutputs.Pdu          - Delta-u constratins
%           secOutputs.Ns           - number of iterations for MRPIS
% OUTPUTS:
%           secOutputs.Tube         - minimal robust positive invariant set
%           secOutputs.Px_robust    - robustustified state constraint
%           secOutputs.Pu_robust    - robustustified input constraint
%           secOutputs.Pdu_robust   - robustustified delta input constraint


% Robust positive invariant set
Tube = ConstructExplicitTube(model,secOutputs);
if Tube.isEmptySet
    error('Cannot continue as the minimal positive invariant set is empty!')
end

% State constraints robustification
Px_robust = secOutputs.Px - Tube;
if Px_robust.isEmptySet
    error('Cannot continue as the state constraint set is empty!')
end
% Input constraints robustification
UKtube = -secOutputs.K*Tube;
Pu_robust = secOutputs.Pu - UKtube;
if Pu_robust.isEmptySet
    error('Cannot continue as the input constraint set is empty!')
end
% Delta input constraints robustification
if any([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
    Pdu_robust = secOutputs.Pdu - UKtube;
end

% Store all computed sets
secOutputs.Tube = Tube;
secOutputs.Px_robust = Px_robust;
secOutputs.Pu_robust = Pu_robust;
if any([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
    secOutputs.Pdu_robust = Pdu_robust;
end
end



function [secOutputs] = DefineMPCProblem_ImplicitTube(secOutputs,model,options)

% Determine MPC formulations that were given by the user
MPCforms = DeterminTEMPCFormulations(model);


% Define decision variables - YALMIP variables
X0 = sdpvar(model.nx,1);
u = sdpvar(model.nu,secOutputs.N,'full');
x = sdpvar(model.nx,secOutputs.N+1+secOutputs.Nz,'full');
w = sdpvar(model.nx,secOutputs.Ns,'full');


% Initialize constraints and objective function
con = [];
obj = 0;


% ------------------------------ Constraints ------------------------------
% Initial condition:      X0 \in X
con = [con, secOutputs.Px.A*X0 <= secOutputs.Px.b ];

% Disturbances:         w \in W
for k = 1:secOutputs.Ns
    con = [con, secOutputs.Pw.A*w(:,k) <= secOutputs.Pw.b];
end

% Implicit tube:
%    ImplicitTube = (1 - alpha)^(-1)*sum_{j=1}^{Ns}((A + B*K)^(j-1)*w(:,j))
%    x0 - x(:,1) == ImplicitTube
SUMA = 0;
for k = 1:length(w)
    SUMA = SUMA + (model.A + model.B*secOutputs.K)^(k-1)*w(:,k);
end
ImplicitTube = (1 - secOutputs.alpha)^(-1)*SUMA;
con = [ con,  X0 == x(:,1) + ImplicitTube];

% Feasibility along the prediction horizon
for k = 1:secOutputs.N
    % Nominal LTI model
    con = [con, x(:,k+1) == model.A*x(:,k) + model.B*u(:,k) ];

    % State-Input constraints
    con = [con, secOutputs.Pxu.A*[x(:,k);u(:,k)] <= secOutputs.Pxu.b - secOutputs.f];
end

% Delta-U constraints
if MPCforms.deltaucon
    uprev = sdpvar(model.nu,1,'full');
    Pdu = secOutputs.Pdu;
    fu = secOutputs.f(end-size(secOutputs.Pu.A,1)+1:end);
    for k = 1:secOutputs.N
        if k == 1
            con = [con, secOutputs.Pu.A*(u(:,k) - uprev) <= secOutputs.Pu.b - fu];
        else
            con = [con, secOutputs.Pu.A*(u(:,k) - u(:,k-1)) <= secOutputs.Pu.b - fu];
        end
    end
end



% ------------------------------ Objective --------------------------------
% We assume following objective function:
%   J = J_terminal + J_Nz + \sum_{k=1}^{N}(x(:,k)'*Qx*x(:,k) + u(:,k)'*R*u(:,k))
% where J_terminal is terminal penalty and J_Nz is "Nz" penalty defined below.
for k = 1:secOutputs.N
    % Objective - state
    obj = obj + x(:,k)'*model.x.penalty.weight*x(:,k);
    % Objective - inputs
    obj = obj + u(:,k)'*model.u.penalty.weight*u(:,k);
end

% --------------------------- Terminal obj-set ----------------------------
% We assume following "Nz" penalty, i.e. where Kz controller is used,
%           J_Nz =  (x(:,N+Nz+1)'*P*x(:,N+Nz+1)) + 
%                   + \sum_{k=N+1}^{Nz}(x(:,k)'*Qz*x(:,k))
% and following "Nz" constraints:
%      x(:,N+k+1) = (A + B*Kz)*x(:,N+k),            k = 1:Nz
%       [x(:,N+k), Kz*x(:,N+k)] \in  XU,            k = 1:Nz 
% where XU are state-input constraints. Finally, the terminal penalty is:
%                   J_terminal = x(:,end)'*Pz*x(:,end)
% and terminal constraint:
%               [x(:,end);Kz*x(:,end)] \in XU - Tube      if         Nz > 0
%               [x(:,end);K*x(:,end)] \in XU - Tube       if         Nz = 0

% LQR stability:
% Introduce LQR-based terminal set and terminal penalty evaluated for LTI
% system, not ULTI system.
if options.LQRstability
    model = ComputeLQRSet(model);
end

if secOutputs.Nz
    Qz = model.x.penalty.weight + secOutputs.Kz'*model.u.penalty.weight*secOutputs.Kz;
    if any(eig(Qz') <= 0), error('MPT+: Ups, penalty matrix Qz is not positive def!'); end

    for k = 1:secOutputs.Nz
        % Complete state predictions
        con = [con, x(:,secOutputs.N+k+1) == (model.A + model.B*secOutputs.Kz)*x(:,secOutputs.N+k)];

        % Penalty
        obj = obj + x(:,secOutputs.N+k)'*Qz*x(:,secOutputs.N+k);

        % Constraints
        con = [con, secOutputs.Pxu.A*[x(:,secOutputs.N+k);secOutputs.Kz*x(:,secOutputs.N+k)] <= secOutputs.Pxu.b - secOutputs.f];
    end
end

% Terminal penalty
if isprop(model.x,'terminalPenalty')
    if secOutputs.Nz
        terminalPenalty = secOutputs.Pz;
    else
        terminalPenalty = model.x.terminalPenalty.weight;
    end
    obj = obj + x(:,end)'*terminalPenalty*x(:,end);
end

% Terminal constraint
if isprop(model.x,'terminalSet')
    if secOutputs.Nz
        con = [con, secOutputs.Pxu.A*[x(:,end);secOutputs.Kz*x(:,end)] <= secOutputs.Pxu.b - secOutputs.f];
    else
        con = [con, secOutputs.Pxu.A*[x(:,end);secOutputs.K*x(:,end)] <= secOutputs.Pxu.b - secOutputs.f];
    end
end

% ---------------------------- MPC definition -----------------------------
% Optimization settings - silent verbose mode
opt_settings = sdpsettings('verbose',0);

p_required = X0;
if any([MPCforms.deltaucon,MPCforms.deltaupen])
    p_required = [p_required; uprev];
end
% end

% Note that the TEMPC will have always p_return = [u(:,1); x(:,1)].
% If required, by solType ==1, then this primal function will be
% overwritten inside of the .toExplicit function
OptimProblem.options = options;
OptimProblem.MPCforms = MPCforms;
OptimProblem.MPCproblem = Opt(con, obj, p_required, ...
    [u(:,1); x(:,1)]);
% construct compact/exact optimizers
OptimProblem.MPCoptimizers.compact = optimizer(con, obj, opt_settings, ...
    p_required, [obj(:); u(:,1) + secOutputs.K*( X0 - x(:,1))]);
OptimProblem.MPCoptimizers.exact = optimizer(con, obj, opt_settings, ...
    p_required, [obj(:); u(:); x(:)]);
% Save all results
secOutputs.OptimProblem = OptimProblem;
end


function [secOutputs] = DefineMPCProblem_ExplicitTube(secOutputs,model,options)
% Function defines optimization problem of an MPC policy
%
% INPUTS:
%     secOutputs.N          - Prediction horizon
%     secOutputs.Px         - original state consraint
%     secOutputs.Tube       - minimal robust invariant set
%     secOutputs.Pw         - disturbance set
%     secOutputs.Px_robust  - robustustified state constraint
%     secOutputs.Pu_robust  - robustustified input constraint
%     secOutputs.K          - LQR gain (u = K*x)
%     model                 - defined model as MPT3 object
%     options.LQRstability  - {0} Compute the terminal penalty/set w.r.t.
%                             the ULTI system (given in 'model')
%                             {1} Compute the terminal penalty/set w.r.t.
%                             the LTI system (not the ULTI as in 'model')
%
% OUTPUTS:
% secOutputs.OptimProblem.MPCproblem                 - optimisation problem (mpt3)
% secOutputs.OptimProblem.Constraints                - opt constraints (YALMIP)
% secOutputs.OptimProblem.Objective                  - opt objective function (YALMIP)
% secOutputs.OptimProblem.Param_Required             - opt required parameters (YALMIP)
% secOutputs.OptimProblem.Param_Returned             - opt returned parameters (YALMIP)
% secOutputs.OptimProblem.CompactForm_Param_Returned - opt required parameters (YALMIP)
% secOutputs.OptimProblem.ExactForm_Param_Returned   - opt returned parameters (YALMIP)
% secOutputs.OptimProblem.YALMIP_settings            - opt settings (YALMIP)
% secOutputs.OptimProblem.options                    - tube MPC options
% secOutputs.OptimProblem.MPCoptimizers.implicit     - optimizer that retuns
%                                                      Utube = u_opt +
%                                                      + K*(x-x_opt) (YALMIP)
% secOutputs.OptimProblem.MPCoptimizers.explicit     - optimizer that retuns
%                                                      [u_opt, x_opt]  (YALMIP)


% Determine MPC formulations that were given by the user
MPCforms = DeterminTEMPCFormulations(model);


% Define decision variables - YALMIP variables
X0 = sdpvar(model.nx,1);
u = sdpvar(model.nu,secOutputs.N,'full');
x = sdpvar(model.nx,secOutputs.N+1+secOutputs.Nz,'full');

% Initialize constraints and objective function
con = [];
obj = 0;

% ------------------------------ Constraints ------------------------------
% We consider following constraints:
%  		 X0 \in X
% 		(X0 - x(:,1)) \in TUBE
%		x(:,k+1) = A*x(:,k) + B*u(:,k)		for k = 1,...,N
%		x(:,k) \in X - TUBE					for k = 1,...,N
%		u(:,k) \in X - K*TUBE			    for k = 1,...,N
%		(u(:,k) - u(:,k-1)) \in dU			for k = 1,...,N
%
%       x(:,N+k+1) = A*x(:,N+k) + B*Kz		for k = 1,...,N
%		x(:,N+k) \in X - TUBE				for k = 1,...,N
%		u(:,N+k) \in X - K*TUBE			    for k = 1,...,N

% Initial condition:      X0 \in X
con = [con, secOutputs.Px.A*X0 <= secOutputs.Px.b ];

% Forcing the "nominal" state x to be inside the tube (centered by X0)
con = [con, secOutputs.Tube.A*( X0 - x(:,1) ) <= secOutputs.Tube.b ];

% Feasibility along the prediction horizon N
for k = 1:secOutputs.N
    % Nominal LTI model
    con = [con, x(:,k+1) == model.A*x(:,k) + model.B*u(:,k) ];

    % Input and state constraints
    con = [con, secOutputs.Px_robust.A*x(:,k) <= secOutputs.Px_robust.b ];
    con = [con, secOutputs.Pu_robust.A*u(:,k) <= secOutputs.Pu_robust.b ];
end

% Delta-U constraints
if MPCforms.deltaucon
    uprev = sdpvar(model.nu,1,'full');
    Pdu = secOutputs.Pdu;
    for k = 1:secOutputs.N
        if k == 1
            con = [con, Pdu.A*(u(:,k) - uprev) <= Pdu.b];
        else
            con = [con, Pdu.A*(u(:,k) - u(:,k-1)) <= Pdu.b];
        end
    end
end

if secOutputs.Nz
    for k = 1:secOutputs.Nz
        % Complete state predictions
        con = [con, x(:,secOutputs.N+k+1) == (model.A + model.B*secOutputs.Kz)*x(:,secOutputs.N+k)];

		% Input and state constraints
		con = [con, secOutputs.Px_robust.A*x(:,secOutputs.N+k) <= secOutputs.Px_robust.b ];
		con = [con, secOutputs.Pu_robust.A*secOutputs.Kz*x(:,secOutputs.N+k) <= secOutputs.Pu_robust.b ];
    end
end

% ------------------------------ Objective --------------------------------
% We assume following objective function:
%   				J = J_terminal + J_Nz + J_nominal 
% where
% 		J_nominal = \sum_{k=1}^{N}    (x(:,k)'*Qx*x(:,k) + u(:,k)'*R*u(:,k))
% 		     J_Nz = \sum_{k=N+1}^{Nz} (x(:,k)'*Qz*x(:,k))
% 	   J_terminal = x(:,end)'*Pz*x(:,end) 		if   Nz > 0
% 	   J_terminal = x(:,end)'*P*x(:,end) 		if   Nz = 0
% with 
% 		Qz = Q + Kz'*R*Kz


% J_nominal
for k = 1:secOutputs.N
    % Objective - state
    obj = obj + x(:,k)'*model.x.penalty.weight*x(:,k);
    % Objective - inputs
    obj = obj + u(:,k)'*model.u.penalty.weight*u(:,k);
end

% J_Nz
if secOutputs.Nz
    Qz = model.x.penalty.weight + secOutputs.Kz'*model.u.penalty.weight*secOutputs.Kz;
    if any(eig(Qz') <= 0), error('MPT+: Ups, penalty matrix Qz is not positive def!'); end

    for k = 1:secOutputs.Nz
        % Penalty
        obj = obj + x(:,secOutputs.N+k)'*Qz*x(:,secOutputs.N+k);
    end
end

% --------------------------- Terminal obj-set ----------------------------
% LQR stability:
% Introduce LQR-based terminal set and terminal penalty evaluated for LTI
% system, not ULTI system.
if options.LQRstability
    model = ComputeLQRSet(model);
end

% Terminal penalty
if isprop(model.x,'terminalPenalty')
    if secOutputs.Nz
        terminalPenalty = secOutputs.Pz;
    else
        terminalPenalty = model.x.terminalPenalty.weight;
    end
    obj = obj + x(:,end)'*terminalPenalty*x(:,end);
end

% Terminal constraint
if isprop(model.x,'terminalSet')
    if secOutputs.Nz
		con = [con, secOutputs.Px_robust.A*x(:,end) <= secOutputs.Px_robust.b ];
		con = [con, secOutputs.Pu_robust.A*secOutputs.Kz*x(:,end) <= secOutputs.Pu_robust.b ];
    else
        con = [con, model.x.terminalSet.A*x(:,end) <= model.x.terminalSet.b ];
    end
end

% ---------------------------- MPC definition -----------------------------
% Optimization settings - silent verbose mode
opt_settings = sdpsettings('verbose',0);

p_required = X0;
if any([MPCforms.deltaucon,MPCforms.deltaupen])
    p_required = [p_required; uprev];
end

% Note that the TEMPC will have always p_return = [u(:,1); x(:,1)].
% If required, by solType ==1, then this primal function will be
% overwritten inside of the .toExplicit function
OptimProblem.options = options;
OptimProblem.MPCforms = MPCforms;
OptimProblem.MPCproblem = Opt(con, obj, p_required, ...
    [u(:,1); x(:,1)]);
% construct compact/exact optimizers
OptimProblem.MPCoptimizers.compact = optimizer(con, obj, opt_settings, ...
    p_required, [obj(:); u(:,1) + secOutputs.K*( X0 - x(:,1))]);
OptimProblem.MPCoptimizers.exact = optimizer(con, obj, opt_settings, ...
    p_required, [obj(:); u(:); x(:)]);
% Save all results
secOutputs.OptimProblem = OptimProblem;
end


function [model] = ComputeLQRSet(model)
Qz = model.x.penalty.weight;

% Import parameters from ULTI system to LTI system
LTI_model = LTISystem('A', model.A, 'B', model.B);
LTI_model.x.penalty = QuadFunction( Qz );
LTI_model.u.penalty = QuadFunction( model.u.penalty.weight );
% Import constraints
if model.x.hasFilter('setConstraint')
    LTI_model.x.setConstraint = model.x.setConstraint;
else
    LTI_model.x.min = model.x.min;
    LTI_model.x.max = model.x.max;
end
if model.u.hasFilter('setConstraint')
    LTI_model.u.setConstraint = model.u.setConstraint;
else
    LTI_model.u.min = model.u.min;
    LTI_model.u.max = model.u.max;
end
% Generate the terminal set
LTI_model.x.with('terminalSet');
LTI_model.x.terminalSet = LTI_model.LQRSet;
% Generate the terminal penalty
[ ~, P ] = dlqr( model.A, model.B, model.x.penalty.weight, model.u.penalty.weight );
LTI_model.x.with('terminalPenalty');
LTI_model.x.terminalPenalty = QuadFunction( P );
% Export parameters from LTI system to ULTI system
% Note, the original termianl set and termianl penalties are
% oberwirtten
if model.x.hasFilter('terminalPenalty')
    fprintf('Warning: The original terminal penalty is overwritten.\n');
else
    model.x.with('terminalPenalty');
end
if model.x.hasFilter('terminalSet')
    fprintf('Warning: The original terminal set is overwritten.\n');
else
    model.x.with('terminalSet');
end
model.x.terminalSet = LTI_model.x.terminalSet;
model.x.terminalPenalty = LTI_model.x.terminalPenalty;
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



function Tube = ConstructExplicitTube(model,secOutputs)
% Invariant Approximations of the Minimal Robust Positively Invariant Set
% by S. V. Rakovic, E. C. Kerrigan, K. I. Kouramas, D. Q. Mayne
% IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 50, NO. 3, MARCH 2005
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1406138

% Extract data
A = model.A;
B = model.B;
K = secOutputs.K;
W = secOutputs.Pw;
Ns = secOutputs.Ns;
alpha = secOutputs.alpha;

% Closed-loop/autonomous system matrix
Acl = A + B*K;

% Eq(2): Fs = Minkowski_Sum_{i=0}^{s-1}( Acl^(i)*W ), where F_0 = {0}
tic;
Fs = W;
for i = 1 : (Ns - 1)
    if size(Fs.V,1) > 2e5
        input('\nMPT+: The tube is too complex! \nPress enter to continue or terminate the script.\n(Consider using "TubeType = implicit" instead of the "TubeType = explicit")')
    end
    ShowProgress('Construction of an explicit tube...',i,Ns-1)
    Fs = Fs + Acl^(i)*W;
    Fs.minVRep();
end
Tube = (1 - alpha)^(-1)*Fs;
compTime = toc;
% fprintf('...done! (computation time: %.2f, #iterations: %.0f, #halfspaces = %.0f)\n',compTime,s,size(Fs.H,1))
fprintf('\b (computation time: %.2f, #iterations: %.0f, #halfspaces = %.0f)\n',compTime,Ns,size(Fs.H,1))
end


function  [form] = DeterminTEMPCFormulations(model)
% Available MPC forms
form.deltaucon = 0;
form.deltaupen = 0;
form.ytrack = 0;
form.xtrack = 0;

if  any([isprop(model.u,'deltaMin'),isprop(model.u,'deltaMax')])
    form.deltaucon = 1;
    fprintf('MPC formulation: Delta-u constraints detected!\n')
end

if isprop(model.u,'deltaPenalty')
    form.deltaupen = 1;
    error('MPTplus does not support Tube MPC with delta-u penalization, yet.')
end

if isprop(model.x,'reference')
    form.xtrack = 1;
    error('MPTplus does not support Tube MPC for state reference tracking problem, yet.')
end

if isprop(model.y,'reference')
    form.ytrack = 1;
    error('MPTplus does not support Tube MPC for output reference tracking problem, yet.')
end
end



function [Ns, alpha, alpha_tol ] = alpha_Ns(model,secOutputs,tol,MaxIter)
% Invariant Approximations of the Minimal Robust Positively Invariant Set
% by S. V. Rakovic, E. C. Kerrigan, K. I. Kouramas, D. Q. Mayne
% IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 50, NO. 3, MARCH 2005
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1406138

% Extract data
A = model.A;
B = model.B;
K = secOutputs.K;
W = secOutputs.Pw;

if (nargin <= 2) 
    tol = 1e-4; % Default tolerance 
end % if 
if ( tol <= 0 )
	tol = 1e-4; % Default tolerance 
end % if 

% Closed-loop/autonomous system matrix 
Acl = A + B*K;

% Problem size
nx = size(Acl, 1); 

% Initial values
Ns = -1; % "s" is integer going toward zero: s -> 0
alpha = 1e6; % Initialize "alpha"
Ms = 1000;
mss = zeros(2*nx,1); % Initialize evaluation of Eq(13)

while( alpha > tol/( tol + Ms ) ) && ( Ns < MaxIter )
    % Increment "s"
    Ns = Ns + 1;
%     fprintf('Iteration of Ns: %d / %d\n',Ns,MaxIter)
    
    % Eq(10): (Acl^s)*W \subseteq alpha*W == h_W( (Acl^(s))'*f_i ) <= alpha*g_i
    % Eq(11): alpha(s) = max_{i \in I} ( ( h_W( (Acl^(s))'*f_i ) ) / g_i )
    alpha = max(W.support(Acl^(Ns)*(W.A)')./W.b);
    
    % Eq(13): M(s) = max_{j=1,...,n} ( sum_{i=0}^{s-1} ( h_W( (Acl^(i))'*e_j ) ), sum_{i=0}^{s-1} ( h_W( -(Acl^(i))'*e_j ) ) )
    mss = mss + W.support([Acl^(Ns), -Acl^(Ns)]);
    Ms = max(mss);
end
alpha_tol = tol/( tol + Ms ); % final tolerance for ALPHA-convergence
end


function f = support_f(model,secOutputs,Ns,alpha)

% Extract inputs:
A = model.A;
B = model.B;
K = secOutputs.K;
W = secOutputs.Pw;
c = secOutputs.Pxu.A(:,1:model.nx);
d = secOutputs.Pxu.A(:,model.nx+1:end);

ni = []; 
f = [];
for i = 1:size(d,1)
    ni = [ni, c(i,:)'+K'*d(i,:)'];
    suma_f = 0;
    for j = 0:Ns-1
        suma_f = suma_f+W.support(((A+B*K)^j)'*ni(:,i));
    end
    f = [f; (1-alpha)^(-1)*suma_f];
end
end
