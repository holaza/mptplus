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
%     TubeType              - {'implicit', 'explicit','elastic',
%                              'timevarying'} Type of the tube embedded in 
%                             the optimization problem. By default, 
%                             TubeType = 'explicit'.
%     Nz                    - number of terminal steps. By default, Nz = 0.
%                             (i.e. prediction steps with control actions 
%                              computed as u_k = Kz*x, where Kz is also 
%                              given by the user)
%     Kz                    - terminal LQR 
%     Pz                    - terminal cost (when Kz is used) 
%     qa                    - weighting constant to penalize the elasticity
%                             of the tube (used only when 
%                             TubeType = 'elastic')
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
if nargin < 2, error('MPTplus: Ups, not enough inputs were selected'); end
% We do not support parametric uncertainties
if iscell(model.A) || iscell(model.B)
    error('MPTplus: Implemented Tube MPC design method does not support parametric uncertainties (LPV systems), yet.')
end

% Define option validation handles
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
isBinary = @(x) ismember(x,[0 1]);
isValidMatrix = @(x) all([all(~isnan(x(:))),all(isnumeric(x(:)))]);
isValidTubeType = @(x) any(strcmpi(x,{'explicit','implicit','elastic','timevarying'}));

% set options/default values
ip = inputParser;
% ip.addParameter('TEMPC', 0);
ip.addParameter('adjList', 0, isBinary);                        % General
ip.addParameter('Tube_tol', 1e-4, validScalarPosNum);           % General
ip.addParameter('Tube_MaxIter', 1e3, validScalarPosNum);        % General
ip.addParameter('solType', 1, isBinary);                        % General
ip.addParameter('LQRstability', 0, isBinary);                   % General
ip.addParameter('TubeType', 'explicit', isValidTubeType);       % General
ip.addParameter('Nz', 0, validScalar);                          % Implicit TMPC 
ip.addParameter('Kz', nan, isValidMatrix);                      % Implicit TMPC / Elastic TMPC
ip.addParameter('Pz', nan, isValidMatrix);                      % Implicit TMPC
ip.addParameter('qa', 1e0, validScalarPosNum);                  % Elastic TMPC
ip.addParameter('forceRigidTubeInETMPC', 0, isBinary);          % Elastic TMPC
if nargin == 3, ip.parse(option{:}); else ip.parse(); end
options = ip.Results;

%% Check if the user selected feasible MPC setup
CheckFeasibleSetup(model,options)

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
    secOutputs = DefineMPCProblem_ImplicitTube(workingModel,secOutputs,options);
elseif strcmpi(options.TubeType,'explicit')
    secOutputs = PreProcessing_ExplicitTube(workingModel,secOutputs);
    secOutputs = DefineMPCProblem_ExplicitTube(workingModel,secOutputs,options);
elseif strcmpi(options.TubeType,'elastic')
    secOutputs = PreProcessing_ElasticTube(workingModel,secOutputs,options);
    secOutputs = DefineMPCProblem_ElasticTube(workingModel,secOutputs,options);
elseif strcmpi(options.TubeType,'timevarying')
%     https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10694700
    secOutputs = PreProcessing_TimeVaryingTube(workingModel,secOutputs);
    secOutputs = DefineMPCProblem_TimeVaryingTube(workingModel,secOutputs,options);
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
if any(P.b <= 0), error('MPTplus: Ups, normalization of the polyhedron failed!'); end
PP = Polyhedron('A',P.A.*P.b.^(-1),'b',P.b.*P.b.^(-1));
end



function [secOutputs] = PreprocessOptions(model,secOutputs,options)
secOutputs.Nz = options.Nz;
secOutputs.Kz = options.Kz;
secOutputs.Pz = options.Pz;

% Is Nz valid
if any([isnan(options.Nz),~isnumeric(options.Nz)])
    error('MPTplus: Parameter "Nz" must be numeric')
elseif options.Nz > secOutputs.N
    error('MPTplus: Parameter "Nz" cannot be larger than the prediciton horizon "N"!')
elseif options.Nz == secOutputs.N
    warning('MPTplus: Parameter "Nz" is equal to the prediction horizon "N".')
end


% Check Kz, Pz inputs only if Nz > 0
if options.Nz > 0
    % Both Kz and Pz need to be given
    if any([sum(isnan(options.Kz(:))),sum(isnan(options.Pz(:)))])
        issuestr = 'MPTplus: not enough inputs were selected. When requesting Nz > 0, then it is also needed to provide: \n';
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
        error('MPTplus: Provided matrix "Kz" has invalid dimensions!')
    end
    [t1,t2] = size(options.Pz);
    if any([t1 ~= t2,t2 ~= model.nx])
        error('MPTplus: Provided matrix "Pz" has invalid dimensions!')
    end

    % Is the Kz stabilizing?
    if any(abs(eig(model.A + model.B*options.Kz)) >= 1)
        warning('MPTplus: Non-stabilizing LQR gain "Kz" was provided, i.e., (abs(eig(A+B*Kz))<1) does not hold.')
    end

    % Is the Pz positive semi-definite?
    if any(eig(options.Pz) < 0)
        warning('MPTplus: Provided terminal penalty "Pz" is not semi-positive define, i.e., all(eig(options.Pz) <= 0) does not hold.')
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
    error('MPTplus: Cannot continue as the state constraint set is empty!')
end
% Input constraints robustification
Pu_robust = Polyhedron('A',secOutputs.Pu.A,'b',secOutputs.Pu.b - fu);
if Pu_robust.isEmptySet
    error('MPTplus: Cannot continue as the input constraint set is empty!')
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



function [secOutputs] = PreProcessing_TimeVaryingTube(model,secOutputs)
% Function preprocess variables for further processing.
%
% INPUTS:
%           model                   - model of the system (as MPT3 object)
%           secOutputs.K            - LQR gain
%           secOutputs.Pw           - Bounded additive uncertainty
%           secOutputs.Px           - state constraints 
%           secOutputs.Pu           - input constraints
%           secOutputs.N            - Prediction horizon
% OUTPUTS:
%           secOutputs.params_TVTMPC.Xf - maximal robust positive invariant set



% Checks of supported functionalities:
if secOutputs.Nz
    error(['MPTplus: We are sorry, but multiple terminal steps Nz' ... 
        'are not supported in the timevarying tube MPC formulation.']);
end
% if any([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
%     error('MPTplus: Ups, delta constraints are not supported!')
% end


% Robust positive invariant set
% Xf = maximal_robust_invariant_set(model,secOutputs);
out = compute_maximal_positive_invariant_set(model, secOutputs);
Xf = out.R_max;
if Xf.isEmptySet
    error('Cannot continue as the minimal positive invariant set is empty!')
end

% Store all computed
secOutputs.params_TVTMPC.Xf = Xf;
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



function [secOutputs] = DefineMPCProblem_ImplicitTube(model,secOutputs,options)

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

% % Delta-U constraints
% if MPCforms.deltaucon
%     uprev = sdpvar(model.nu,1,'full');
%     Pdu = secOutputs.Pdu;
%     fu = secOutputs.f(end-size(Pdu.A,1)+1:end);
%     for k = 1:secOutputs.N
%         if k == 1
%             con = [con, Pdu.A*(u(:,k) - uprev) <= Pdu.b - fu];
%         else
%             con = [con, Pdu.A*(u(:,k) - u(:,k-1)) <= Pdu.b - fu];
%         end
%     end
% end
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
    if any(eig(Qz') <= 0), error('MPTplus: Ups, penalty matrix Qz is not positive def!'); end

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



function [secOutputs] = DefineMPCProblem_ExplicitTube(model,secOutputs,options)
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
    if any(eig(Qz') <= 0), error('MPTplus: Ups, penalty matrix Qz is not positive def!'); end

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
        input('\nMPTplus: The tube is too complex! \nPress enter to continue or terminate the script.\n(Consider using "TubeType = implicit" instead of the "TubeType = explicit")')
    end
    ShowProgress('Construction of an explicit tube...',i,Ns-1)
    Fs = Fs + Acl^(i)*W;
    Fs.minVRep();
end
Tube = (1 - alpha)^(-1)*Fs;
compTime = toc;
% fprintf('...done! (computation time: %.2f, #iterations: %.0f, #halfspaces = %.0f)\n',compTime,s,size(Fs.H,1))
fprintf('\b (computation time: %.2f, #iterations: %.0f, #halfspaces = %.0f)\n',compTime,Ns,size(Fs.H,1))

% Normalize the tube:
Tube = NormalizePolyhedron(Tube.minHRep);
end


function  [form] = DeterminTEMPCFormulations(model)
% Available MPC forms
form.deltaucon = 0;
form.deltaupen = 0;
form.ytrack = 0;
form.xtrack = 0;
form.xterminalPenalty = 0;
form.ublock = 0;

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

if isprop(model.x,'terminalPenalty')
    form.xterminalPenalty = 1;
end

if isprop(model.u,'block')
    form.ublock = 1;
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



function [secOutputs] = PreProcessing_ElasticTube(model,secOutputs,options)
% -------------------------------------------------------------------------
%       [secOutputs] = PreProcessing_ElasticTube(model,secOutputs,options)
% -------------------------------------------------------------------------
% This function computes approximation parameters 
% L - Global Elastic Tubes Dynamics (C*A*z + C*B*v + phi(sigma) + d <= C*z+ + a+)
%     approximated by C*A*z + C*B*v + L*sigma + d <= C*z+ + a+
% M - Global Elastic Tubes State Constraints (G*z + gamma(sigma) <= 1) 
%     approximated by (G*z + M*sigma <= 1)
% T - Global Elastic Tubes Control Constraints (H*v + ni(sigma) <= 1)
%     approximated by (H*v + T*sigma <= 1)
% where z are nominal states and v are nominal control inputs.      
% -------------------------------------------------------------------------
% INPUTS:
%       model.A                 - system dynamics: x+ = A*x + B*u
%       model.B                 - system dynamics: x+ = A*x + B*u
%       model.x.penalty.weight  - state penalization
%       model.u.penalty.weight  - input penalization
%       secOutputs.K            - LQR policy
%       secOutputs.Px           - state constraints (normalized), i.e.
%                                      X := {x \in R^n | G*x <= 1}
%       secOutputs.Pu           - input constraints (normalized), i.e.
%                                      U := {u \in R^m | H*u <= 1}
%       secOutputs.Pw           - disturbance constraints (normalized),i.e.
%                                      W := {w \in R^n | D*w <= 1}
%       options.qa              - weighting constant for sigma penalization
%                                 (sigma is the optimized variable that 
%                                  will modify the tube C*x <= 1 as: 
%                                  C*x <= sigma)
%       options.Kz              - terminal LQR gain
%
% OUTPUTS:
%       L(1) \in R^{qs \times qs}_{>=0}, such that L*C = C*(A + B*K(sigma))
%       M(1) \in R^{qx \times qs}_{>=0}, such that M*C = G
%       T(1) \in R^{qu \times qs}_{>=0}, such that T*C = H*K(sigma)
%       d = {d_i \in R | d_i = support(W,C_i), i = 1, ..., qs}
%       gamma(1): R^{qs}_{>=0} -> R^{qx}_{>=0}
%       ni(1): R^{qs}_{>=0} -> R^{qu}_{>=0}
%       fi(1): R^{qs}_{>=0} -> R^{qs}_{>=0}
%       Tset - terminal constraint
%       Px_terminal - terminal state penalty
%       Pa_terminal - terminal sigma penalty
%       Qa          - penalization of the sigma variables
%
% Reference:
% Sasa V.Rakovic, William S. Levine and Behcet Acikmese, Elastic Tube Model
% Predictive Control, ACC 2016.

% Notation Background:
%       Inputs:
%       system dynamics matrices model.A,model.B as: x+ = A*x + B*u
%       X := {x \in R^n | G*x <= 1}
%       U := {u \in R^m | H*u <= 1}
%       W := {w \in R^n | D*w <= 1}
%       S(sigma) := {sigma \in R^n | C*x <= sigma}
%       K(sigma): R^{qs} -> R^{n \times m}
%       Outputs:
%       L \in R^{qs \times qs}_{>=0}, such that L*C = C*(A + B*K(sigma))
%       M \in R^{qx \times qs}_{>=0}, such that M*C = G
%       T \in R^{qu \times qs}_{>=0}, such that T*C = H*K(sigma)
%       d = {d_i \in R | d_i = support(W,C_i), i = 1, ..., qs}
%       gamma(1): R^{qs}_{>=0} -> R^{qx}_{>=0}
%       ni(1): R^{qs}_{>=0} -> R^{qu}_{>=0}


% To simplify notation
X = secOutputs.Px;
U = secOutputs.Pu;
W = secOutputs.Pw;



% -------------------------------------------------------------------------
% Tube computation (Robust positive invariant set), 
%               S(sigma) := {sigma \in R^n | C*x <= sigma}
% -------------------------------------------------------------------------
Sa = ConstructExplicitTube(model,secOutputs);
if Sa.isEmptySet
    error('Cannot continue as the minimal positive invariant set is empty!')
end



% -------------------------------------------------------------------------
% Dimemnsions:
% -------------------------------------------------------------------------
% qs -> number of cons for the tube
% qx -> number of cons for the state constraints X
% qu -> number of cons for the input constraints U
qs = size(Sa.A,1);
qx = size(X.A,1);
qu = size(U.A,1);

% TODO: In this script we consider sigma to be static and equal to
sigma = ones(qs,1);
% TODO: make this as an optional parameter + make K funtion of the sigma
%               K(sigma): R^{qs} -> R^{n \times m}
Ka = @(sigma) secOutputs.K; 
% TODO: assumption Kz == Klqr == K(1)
if sum(sum(isnan(options.Kz)))
    fprintf('MPTplus: terminal controller not provided, assuming LQR gain as the terminal controller ...\n')
    Kz = Ka(sigma);
else
    Kz = options.Kz;
end
 
% TODO: Possibly make a flaxible penalization of each hyperplane?
Qa = eye(size(Sa.A,1)) * options.qa; 


% Evaluation of inputs:
if options.Nz, warning('MPTplus: Ups, parameters Nz/Pz are not supported in elastic tube MPC. These parameters are omitted ... \n'); end


% -------------------------------------------------------------------------
% Global Elastic Tubes Dynamics: d (3.4)
% -------------------------------------------------------------------------
% Computation of d \in R^{qs}_{>=0} that represents "The biggest impact of 
% disturbances from W to states", i.e.:
% The tube is defiend as Z := {x \in R^n | C*x <= 1}. The disturbance 
% w \in W := {w \in R^n | D*w <= 1} is directly injected to the state: 
% x+ = A*x + B*u + w. From the Tube point of view this distrubance w will 
% have following impact: C*(x + w) <= 1, hence C*x + C*w <= 1
% 
%       d_i = support(W,C_i)       i = 1, ..., qs
%
% MPT3: support(x) := max x'*y s.t. y in Set
d = nan(qs,1);
for k = 1:qs
    d(k) = W.support(Sa.A(k,:)');
end
if any(isnan(d))
    error('MPTplus: Ups, computation of vector d failed!')
end



% -------------------------------------------------------------------------
% Global Elastic Tubes State Constraints:gamma (3.12)
% -------------------------------------------------------------------------
% gamma(sigma): R^{qs}_{>=0} -> R^{qx}_{>=0}

gamma = nan(qx,1);
for k = 1:qx
    gamma(k) = Sa.support(X.A(k,:)');
end
if any(isnan(gamma))
    error('MPTplus: Ups, computation of vector gamma failed!') 
end



% -------------------------------------------------------------------------
% Global Elastic Tubes Control Constraints: ni (3.20)
% -------------------------------------------------------------------------
% n(sigma): R^{qs}_{>=0} -> R^{qu}_{>=0}

KaSa = Ka(sigma)*Polyhedron('A', Sa.A, 'b', sigma);
ni = nan(qu,1);
for k = 1:qu
    ni(k) = KaSa.support(U.A(k,:)');
end
if any(isnan(ni))
    error('MPTplus: Ups, computation of vector ni failed!')
end



% -------------------------------------------------------------------------
% M(1) computation
% -------------------------------------------------------------------------
% M \in R^{qx \times qs}_{>=0}

obj = 0;
con = [];
M = sdpvar(qx,qs,'full');

% Objective function
% min sum_{i=1}^{qx} sum_j^{qs} M(i,j)^2  
for k = 1:qx
    for kk = 1:qs
        obj = obj + M(k,kk)*M(k,kk);
    end
end

% Constraints
% M*1 = gamma(1)
% M >= 0
% M*C = G
con = [con; M*ones(qs,1) == gamma];
con = [con; M >= 0];
con = [con; M*Sa.A == X.A];

opt_settings = sdpsettings('verbose',0);

DIAGNOSTIC = optimize(con,obj,opt_settings);
if DIAGNOSTIC.problem ~= 0
    error('MPTplus: Ups, cant calculate M(1) component!');
end
M_opt = double(M);



% -------------------------------------------------------------------------
% Computation of L(1)
% -------------------------------------------------------------------------
% L \in R^{qs \times qs}_{>=0}

obj = 0;
con = [];
L = sdpvar(qs,qs,'full');

% Objective function
% min sum_{i=1}^{qs} sum_j^{qs} L(i,j)^2  
for k = 1:qs
    for kk = 1:qs
        obj = obj + L(k,kk)*L(k,kk);
    end
end

% Constraints
con = [con; L*ones(qs,1) == ones(qs,1) - d];
con = [con; L >= 0];
con = [con; L*Sa.A == Sa.A*(model.A + model.B*Ka(sigma))];

opt_settings = sdpsettings('verbose',0);

DIAGNOSTIC = optimize(con,obj,opt_settings);
if DIAGNOSTIC.problem ~= 0
    error('MPTplus: Ups, cant calculate L(1) component!');
end
L_opt = double(L);



% -------------------------------------------------------------------------
% Computation of T(1)
% -------------------------------------------------------------------------
% T \in R^{qu \times qs}_{>=0}

obj = 0;
con = [];
T = sdpvar(qu,qs,'full');

% Objective function
% min sum_{i=1}^{qu} sum_j^{qs} T(i,j)^2  
for k = 1:qu
    for kk = 1:qs
        obj = obj + T(k,kk)*T(k,kk);
    end
end

% Constraints
con = [con; T*ones(qs,1) == ni];
con = [con; T >= 0];
con = [con; T*Sa.A == U.A*Ka(sigma)];

opt_settings = sdpsettings('verbose',0);

DIAGNOSTIC = optimize(con,obj,opt_settings);
if DIAGNOSTIC.problem ~= 0
    error('MPTplus: Ups, cant calculate T(1) component!');
end
T_opt = double(T);



% -------------------------------------------------------------------------
% Terminal set Tset
% -------------------------------------------------------------------------
% Tset := {(z,a) \in R^{n+qs} | -a <= 0, G*z + M1*a <= 1, H*Kz*z + T1*a <= 1}

%               -a <= 0
%       G*z + M1*a <= 1
%    H*Kz*z + T1*a <= 1
% ---------------------------
% H*[z a] <= h
H = [zeros(qs,size(model.A,1)) -eye(qs);
    X.A M_opt;
    U.A*Kz T_opt];
h = [zeros(qs,1);ones(qx,1);ones(qu,1)];

% Setup dynamics:
% z+ = (A + B*Kz)*z
% a+ = L1*a + d
% ------------------
% [z+; a+] = [(A + B*Kz) 0; 0 L1] [z; a] + [0;d]


nx = size(model.A,2);
AA = [(model.A + model.B*Kz) zeros(nx,qs); zeros(qs,nx) L_opt];
BB = [zeros(nx,1); d];
model_inv_set = LTISystem('A', AA, 'f', BB);
model_inv_set.x.with('setConstraint');
model_inv_set.x.setConstraint = Polyhedron('A',H, 'b', h);
Tset = model_inv_set.invariantSet();



% -------------------------------------------------------------------------
% Terminal penalty
% -------------------------------------------------------------------------
if all(eig(L_opt) < 1)
    if all(eig(model.A + model.B*Kz) < 1)
        obj = 0;
        con = [];
        Pz = sdpvar(nx);
        Pa = sdpvar(qs);

        Qx = model.x.penalty.weight;
        Qu = model.u.penalty.weight;

        con = [con, (model.A + model.B*Kz)'*Pz*(model.A + model.B*Kz) - ...
            Pz <= -(Qx + Kz'*Qu*Kz)];
        con = [con, L_opt'*Pa*L_opt - Pa <= -Qa];

        opt_settings = sdpsettings('verbose',0);
        DIAGNOSTIC = optimize(con,obj,opt_settings);
        if DIAGNOSTIC.problem ~= 0
            error('MPTplus: Ups, cant calculate terminal penalties Pz and Pa!');
        end
        Px_terminal = double(Pz);
        Pa_terminal = double(Pa);
    else
        error('MPTplus: Ups, unstabilizing terminal controller Kz was selected!')
    end
else
    error('MPTplus: Ups, computed L(1) matrix is wrong!')
end



% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
elasticparams.L1 = L_opt;
elasticparams.M1 = M_opt;
elasticparams.T1 = T_opt;
elasticparams.d = d;
elasticparams.gamma1 = gamma;
elasticparams.ni1 = ni;
elasticparams.fi1 = 1 - d;
elasticparams.qs = qs;
elasticparams.qx = qx;
elasticparams.qu = qu;
elasticparams.sigma = sigma;
elasticparams.Tset = Tset;
elasticparams.Px_terminal = Px_terminal;
elasticparams.Pa_terminal = Pa_terminal;
elasticparams.Qa = Qa;
elasticparams.Ka = Ka;
secOutputs.Tube = Sa;
secOutputs.Px_robust = {};
secOutputs.Pu_robust = {};
secOutputs.elasticparams = elasticparams;

end



function [secOutputs] = DefineMPCProblem_ElasticTube(model,secOutputs,options)
% Function defines optimization problem of an MPC policy
%
% INPUTS:
%     secOutputs.N          - Prediction horizon
%     secOutputs.Px         - original state consraint
%     secOutputs.Tube       - minimal robust invariant set
%     secOutputs.Px_robust  - robustustified state constraint
%     secOutputs.Pu_robust  - robustustified input constraint
%     secOutputs.Pdu        - delta-u constraints
%     secOutputs.K          - LQR gain (u = K*x)
%     model                 - defined model as MPT3 object
%     options.LQRstability  - {0} Compute the terminal penalty/set w.r.t.
%                             the ULTI system (given in 'model')
%                             {1} Compute the terminal penalty/set w.r.t.
%                             the LTI system (not the ULTI as in 'model')
%     options.forceRigidTubeInETMPC  - forces rigit tube (i.e. not elastic)
%     secOutputs.elasticparams - all precomputed elastic parameters (see 
%                                ComputeElasticParams)
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


% TODO:
% !!! in the RHC is used "secOutputs.K" (not secOutputs.elasticparams.Ka)
% !!! Kz/Pz is not impleemnted 



% Determine MPC formulations that were given by the user
MPCforms = DeterminTEMPCFormulations(model);


% Define decision variables - YALMIP variables
X0 = sdpvar(model.nx,1);
u = sdpvar(model.nu,secOutputs.N,'full');
x = sdpvar(model.nx,secOutputs.N+1+secOutputs.Nz,'full');
sigma = sdpvar(secOutputs.elasticparams.qs,secOutputs.N+1+secOutputs.Nz,'full');


% Initialize constraints and objective function
con = [];
obj = 0;

% ------------------------------ Constraints ------------------------------
% We consider following constraints:
%  		 X0 \in X                       <- this is not needed (only imposes infeasibility)
% 		(X0 - x(:,1)) \in S(sigma(:,1))
%       
% 		-sigma(:,k) <= 0                                for k = 1, ..., N
%       C*A*z(:,k) + C*B*v(:,k) + L(1)*sigma(:,k) + 
%       + d - C*z(:,k+1) - sigma(:,k+1) <= 0            for k = 1, ..., N
%       G*z(:,k) + M(1)*sigma(:,k) - 1 <= 0             for k = 1, ..., N
%       H*v(:,k) + T(1)*sigma(:,k) - 1 <= 0             for k = 1, ..., N
%       Fz*z(:,N+1) + Fa*sigma(:,N+1) - Fa*1 - 1 <= 0
%
%       TODO: implement constraints also for Kz > 0

C = secOutputs.Tube.A; % we consider that the tube is normalized 

con = [con, C*X0 - C*x(:,1) - sigma(:,1) <= 0];
for k = 1:secOutputs.N
    con = [con, -sigma(:,k) <= 0];
    % Global elastic tube dynamics
    con = [con, C*model.A*x(:,k) + C*model.B*u(:,k) + secOutputs.elasticparams.L1*sigma(:,k) + ...
           + secOutputs.elasticparams.d - C*x(:,k+1) - sigma(:,k+1) <= 0];
    % Global elasetic tube state constraints
    con = [con, secOutputs.Px.A*x(:,k) + secOutputs.elasticparams.M1*sigma(:,k) - 1 <= 0];
    % Global elasetic tube input constraints
    con = [con, secOutputs.Pu.A*u(:,k) + secOutputs.elasticparams.T1*sigma(:,k) - 1 <= 0];
end

% % Delta-U constraints
% if MPCforms.deltaucon
%     uprev = sdpvar(model.nu,1,'full');
%     Pdu = secOutputs.Pdu;
%     for k = 1:secOutputs.N
%         if k == 1
%             con = [con, Pdu.A*(u(:,k) - uprev) + secOutputs.elasticparams.T1*sigma(:,k) - 1 <= 0];
%         else
%             con = [con, Pdu.A*(u(:,k) - u(:,k-1)) + secOutputs.elasticparams.T1*sigma(:,k) - 1 <= 0];
%         end
%     end
% end
% Delta-U constraints
if MPCforms.deltaucon
    uprev = sdpvar(model.nu,1,'full');
    Pdu = secOutputs.Pdu;
    for k = 1:secOutputs.N
        if k == 1
            con = [con, Pdu.A*(u(:,k) - uprev)  <= Pdu.b];
        else
            con = [con, Pdu.A*(u(:,k) - u(:,k-1)) <= Pdu.b];
        end
    end
end

% TODO: remove this option
if options.forceRigidTubeInETMPC
    for k = 1:secOutputs.N+1
        con = [con, sigma(:,k) == ones(size(sigma(:,k)))];
    end
end


% ------------------------------- Objective -------------------------------
% Consider a stage cost:
%       el(z,v,sigma) = z'*Qz*z + v'*Qv*v + (sigma - 1)'*Qa*(sigma - 1)
% and terminal cost:
%       Vf(z,sigma) = z'*Pz*z + (sigma - 1)'*Pa*(sigma - 1)
% then the overall objective function is given by:
% J(z,v,sigma) = Vf(z(:,N+1),sigma(:,N+1)) + \sum_{k=1}^{N} el(z(:,k),v(:,k),sigma(:,k))

Q = model.x.penalty.weight;
R = model.u.penalty.weight;
Qa = secOutputs.elasticparams.Qa;
for k = 1:secOutputs.N
    obj = obj + x(:,k)'*Q*x(:,k) + u(:,k)'*R*u(:,k) + (1 - sigma(:,k))'*Qa*(1 - sigma(:,k));
end



% --------------------------- Terminal obj/cons ---------------------------
% TODO: let user to force his-her terminal parameters

if options.LQRstability
    model = ComputeLQRSet(model);
end

% Terminal constraints
if isprop(model.x,'terminalSet')
    A = secOutputs.elasticparams.Tset.A;
    b = secOutputs.elasticparams.Tset.b;
    con = [con, A*[x(:,end); sigma(:,end)] <= b];
end

% Terminal penalty
if isprop(model.x,'terminalPenalty')
    Pt = secOutputs.elasticparams.Px_terminal;
    Pa = secOutputs.elasticparams.Pa_terminal;
    N = secOutputs.N;
    obj = obj + x(:,N+1)'*Pt*x(:,N+1) + (1 - sigma(:,N+1))'*Pa*(1 - sigma(:,N+1));
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
    [u(:,1); x(:,1); sigma(:,1)]);
% construct compact/exact optimizers
% OptimProblem.MPCoptimizers.compact = optimizer(con, obj, opt_settings, ...
%     p_required, [obj(:); u(:,1) + secOutputs.K*( X0 - x(:,1))]);
OptimProblem.MPCoptimizers.compact = optimizer(con, obj, opt_settings, ...
    p_required, [obj(:); u(:,1) + secOutputs.elasticparams.Ka(sigma(:,1))*( X0 - x(:,1))]);
% OptimProblem.MPCoptimizers.exact = optimizer(con, obj, opt_settings, ...
%     p_required, [obj(:); u(:); x(:)]);
OptimProblem.MPCoptimizers.exact = optimizer(con, obj, opt_settings, ...
    p_required, [obj(:); u(:); x(:); sigma(:)]);
% Save all results
secOutputs.OptimProblem = OptimProblem;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Supporting functions for  PreProcessing_TimeVaryingTube %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = compute_maximal_positive_invariant_set(model, secOutputs, opts)
% compute_invariant_sets
%
% System:
%   x+ = A x + B u + w
%   u  = K x
%
% Returns:
%   out.LQR_region_max   - maximal nominal invariant set for x+ = (A+B*K)x
%                          subject to x in X and Kx in U
%   out.R_max            - maximal robust positively invariant set for
%                          x+ = (A+B*K)x + w, w in W
%                          subject to x in X and Kx in U
%   out.R_min_outer      - epsilon-outer approximation of minimal RPI set
%   out.X_admissible_cl  - X ∩ {x | Kx ∈ U}
%
% Inputs:
%   X, U, W must be MPT3 Polyhedron objects.
%
% Notes:
%   - This implementation intentionally avoids minHRep(), because that can
%     trigger LP-based simplification routines in MPT3 and fail depending
%     on the LP backend.
%   - The minimal RPI set is returned as an epsilon-outer approximation.
%
% Required:
%   - MPT3
%   - YALMIP
%
% Optional:
%   - LP solver for YALMIP verification

    if nargin < 3
        opts = struct();
    end
    opts = set_default_opts(opts);

    % Extract inputs
    A = model.A;
    B = model.B;
    K = secOutputs.K;
    U = secOutputs.Pu;
    X = secOutputs.Px;
    W = secOutputs.Pw;

    % ------------------------------------------------------------
    % Dimensions and basic checks
    % ------------------------------------------------------------
    nx = size(A,1);
    nu = size(B,2);

    assert(size(A,2) == nx, 'A must be square.');
    assert(size(B,1) == nx, 'A and B dimensions are inconsistent.');
    assert(all(size(K) == [nu, nx]), 'K must be of size nu-by-nx.');

    assert(isa(X, 'Polyhedron'), 'X must be a Polyhedron.');
    assert(isa(U, 'Polyhedron'), 'U must be a Polyhedron.');
    assert(isa(W, 'Polyhedron'), 'W must be a Polyhedron.');

    assert(X.Dim == nx, 'X has wrong dimension.');
    assert(U.Dim == nu, 'U has wrong dimension.');
    assert(W.Dim == nx, 'W has wrong dimension.');

    assert(~X.isEmptySet(), 'X is empty.');
    assert(~U.isEmptySet(), 'U is empty.');
    assert(~W.isEmptySet(), 'W is empty.');

    assert(X.isBounded(), 'X must be bounded.');
    assert(U.isBounded(), 'U must be bounded.');
    assert(W.isBounded(), 'W must be bounded.');

    if ~W.contains(zeros(nx,1))
        error('This implementation assumes 0 ∈ W.');
    end

    Acl = A + B*K;

    % ------------------------------------------------------------
    % Admissible closed-loop state set:
    %   Xad = {x | x ∈ X, Kx ∈ U}
    % ------------------------------------------------------------
    U_pre_K = poly_preimage(U, K);
    Xad = X & U_pre_K;

    if Xad.isEmptySet()
        warning('The admissible closed-loop set X ∩ {x | Kx ∈ U} is empty.');
    end

    % ------------------------------------------------------------
    % 2) Maximal robust positively invariant set under fixed K
    % ------------------------------------------------------------
    [R_max, max_rpi_info] = maximal_rpi_fixed_feedback(Acl, Xad, W, opts);

    % ------------------------------------------------------------
    % Optional verification with YALMIP
    % ------------------------------------------------------------
    verification = struct();

    if opts.verifyWithYALMIP
        verification.R_max_invariant = verify_rpi_yalmip(R_max, Acl, W, opts);
        verification.R_max_admissible = poly_contains_poly(Xad, R_max, opts.tol);
    end

    out = struct();
    out.R_max            = R_max;
    out.X_admissible_cl  = Xad;

    out.info = struct();
    out.info.Acl = Acl;
    out.info.max_rpi = max_rpi_info;
    out.info.verification = verification;
end

% =========================================================================
% OPTIONS
% =========================================================================
function opts = set_default_opts(opts)
    if ~isfield(opts, 'maxIterLQR'),       opts.maxIterLQR = 500; end
    if ~isfield(opts, 'maxIterMaxRPI'),    opts.maxIterMaxRPI = 500; end
    if ~isfield(opts, 'maxIterMinRPI'),    opts.maxIterMinRPI = 1000; end
    if ~isfield(opts, 'epsMRPI'),          opts.epsMRPI = 1e-3; end
    if ~isfield(opts, 'tol'),              opts.tol = 1e-8; end
    if ~isfield(opts, 'verifyWithYALMIP'), opts.verifyWithYALMIP = true; end
    if ~isfield(opts, 'solver'),           opts.solver = ''; end

    assert(opts.maxIterLQR >= 1, 'opts.maxIterLQR must be >= 1.');
    assert(opts.maxIterMaxRPI >= 1, 'opts.maxIterMaxRPI must be >= 1.');
    assert(opts.maxIterMinRPI >= 1, 'opts.maxIterMinRPI must be >= 1.');
    assert(opts.epsMRPI > 0 && opts.epsMRPI < 1, 'opts.epsMRPI must be in (0,1).');
end

% =========================================================================
% PREIMAGE UNDER LINEAR MAP
%   Ppre = {x | Mx ∈ P}
% If P = {y | A y <= b, Ae y = be}, then
%        Ppre = {x | A M x <= b, Ae M x = be}
% =========================================================================
function Ppre = poly_preimage(P, M)
    if isempty(P.Ae)
        Ppre = Polyhedron('A', P.A*M, 'b', P.b);
    else
        Ppre = Polyhedron('A', P.A*M, 'b', P.b, ...
                          'Ae', P.Ae*M, 'be', P.be);
    end
end



% =========================================================================
% 2) MAXIMAL RPI SET FOR x+ = Acl x + w, w ∈ W
%
% Iteration:
%   R_0 = Xad
%   R_{k+1} = Xad ∩ Pre(R_k)
%   Pre(S) = {x | Acl x + w ∈ S, ∀w∈W}
%          = {x | Acl x ∈ S ⊖ W}
% =========================================================================
function [R, info] = maximal_rpi_fixed_feedback(Acl, Xad, W, opts)
    R = Xad;
    converged = false;
    iter = 0;

    for k = 1:opts.maxIterMaxRPI
        iter = k;

        tightened = R - W;   % Pontryagin difference
        if tightened.isEmptySet()
            R = Polyhedron.emptySet(size(Acl,2));
            converged = true;
            break;
        end

        PreR = poly_preimage(tightened, Acl);
        Rnext = Xad & PreR;

        if Rnext.isEmptySet()
            R = Rnext;
            converged = true;
            break;
        end

        if poly_equal(Rnext, R, opts.tol)
            R = Rnext;
            converged = true;
            break;
        end

        R = Rnext;
    end

    info = struct();
    info.iterations = iter;
    info.converged = converged;
end



% =========================================================================
% POLYHEDRON EQUALITY TEST BY MUTUAL CONTAINMENT
% =========================================================================
function tf = poly_equal(P, Q, tol)
    tf = poly_contains_poly(P, Q, tol) && poly_contains_poly(Q, P, tol);
end

% =========================================================================
% CHECK Q ⊆ P
%
% This version prefers a vertex test for bounded Q.
% If vertex enumeration fails, it falls back to YALMIP support-function LPs.
% =========================================================================
function tf = poly_contains_poly(P, Q, tol)
    if Q.isEmptySet()
        tf = true;
        return;
    end

    try
        if Q.isBounded()
            V = Q.V;
            tf = true;
            for i = 1:size(V,1)
                if ~P.contains(V(i,:).', tol)
                    tf = false;
                    return;
                end
            end
            return;
        end
    catch
        % fall through to LP-based test below
    end

    tf = poly_contains_poly_lp(P, Q, tol);
end

% =========================================================================
% LP-BASED CONTAINMENT TEST:
%   Q ⊆ P  iff  max_{x∈Q} a_i x - b_i <= 0 for each facet a_i x <= b_i of P
% =========================================================================
function tf = poly_contains_poly_lp(P, Q, tol)
    if Q.isEmptySet()
        tf = true;
        return;
    end

    H = P.A;
    h = P.b;

    x = sdpvar(Q.Dim,1);
    C = poly_to_yalmip(Q, x, tol);

    ops = sdpsettings('verbose', 0);

    for i = 1:size(H,1)
        obj = H(i,:)*x - h(i);
        sol = optimize(C, -obj, ops);
        if sol.problem ~= 0
            error('Containment LP failed in poly_contains_poly_lp. problem=%d', sol.problem);
        end
        if value(obj) > tol
            tf = false;
            return;
        end
    end

    tf = true;
end

% =========================================================================
% CONVERT POLYHEDRON TO YALMIP CONSTRAINTS
% =========================================================================
function C = poly_to_yalmip(P, z, tol)
    C = [];
    if ~isempty(P.A)
        C = [C, P.A*z <= P.b + tol];
    end
    if ~isempty(P.Ae)
        C = [C, P.Ae*z == P.be];
    end
end

% =========================================================================
% VERIFY NOMINAL INVARIANCE OF P FOR x+ = Acl x
% =========================================================================
function rep = verify_nominal_invariance_yalmip(P, Acl, opts)
    rep = struct('ok', false, 'maxViolation', inf, 'facetValues', []);

    if P.isEmptySet()
        rep.ok = true;
        rep.maxViolation = -inf;
        rep.facetValues = [];
        return;
    end

    H = P.A;
    h = P.b;

    x = sdpvar(size(Acl,2),1);
    C = poly_to_yalmip(P, x, opts.tol);

    ops = sdpsettings('verbose', 0);
    if ~isempty(opts.solver)
        ops = sdpsettings(ops, 'solver', opts.solver);
    end

    vals = nan(size(H,1),1);
    maxViol = -inf;

    for i = 1:size(H,1)
        obj = H(i,:)*(Acl*x) - h(i);
        sol = optimize(C, -obj, ops);
        if sol.problem ~= 0
            error('YALMIP failed during nominal invariance verification. problem=%d', sol.problem);
        end
        vals(i) = value(obj);
        maxViol = max(maxViol, vals(i));
    end

    rep.ok = (maxViol <= opts.tol);
    rep.maxViolation = maxViol;
    rep.facetValues = vals;
end

% =========================================================================
% VERIFY ROBUST POSITIVE INVARIANCE OF P FOR x+ = Acl x + w, w ∈ W
% =========================================================================
function rep = verify_rpi_yalmip(P, Acl, W, opts)
    rep = struct('ok', false, 'maxViolation', inf, 'facetValues', []);

    if P.isEmptySet()
        rep.ok = true;
        rep.maxViolation = -inf;
        rep.facetValues = [];
        return;
    end

    H = P.A;
    h = P.b;

    x = sdpvar(size(Acl,2),1);
    w = sdpvar(size(Acl,1),1);

    C = [];
    C = [C, poly_to_yalmip(P, x, opts.tol)];
    C = [C, poly_to_yalmip(W, w, opts.tol)];

    ops = sdpsettings('verbose', 0);
    if ~isempty(opts.solver)
        ops = sdpsettings(ops, 'solver', opts.solver);
    end

    vals = nan(size(H,1),1);
    maxViol = -inf;

    for i = 1:size(H,1)
        obj = H(i,:)*(Acl*x + w) - h(i);
        sol = optimize(C, -obj, ops);
        if sol.problem ~= 0
            error('YALMIP failed during robust invariance verification. problem=%d', sol.problem);
        end
        vals(i) = value(obj);
        maxViol = max(maxViol, vals(i));
    end

    rep.ok = (maxViol <= opts.tol);
    rep.maxViolation = maxViol;
    rep.facetValues = vals;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [secOutputs] = DefineMPCProblem_TimeVaryingTube(model,secOutputs,options)
% Function defines optimization problem of an MPC policy
%
% INPUTS:
%     secOutputs.N             - Prediction horizon
%     secOutputs.params_TVTMPC - precomputed parameters for the approach
%     secOutputs.Pxu           - combined state-input constraints
%     secOutputs.K             - LQR gain (u = K*x)
%     model                    - defined model as MPT3 object
%     options.LQRstability     - {0} Compute the terminal penalty/set w.r.t.
%                                the ULTI system (given in 'model')
%                                {1} Compute the terminal penalty/set w.r.t.
%                                the LTI system (not the ULTI as in 'model')   
%
% OUTPUTS:
% secOutputs.OptimProblem.MPCproblem                 - optimisation problem (mpt3)
% secOutputs.OptimProblem.options                    - tube MPC options
% secOutputs.OptimProblem.paramtext                  - list of required
%                                                      parameters for the 
%                                                      optimisation problem
% secOutputs.OptimProblem.MPCoptimizers.compact      - optimizer that retuns
%                                                      Utube = u_opt +
%                                                      + K*(x-x_opt) (YALMIP)
% secOutputs.OptimProblem.MPCoptimizers.exact        - optimizer that retuns



% Determine MPC formulations that were given by the user
MPCforms = DeterminTEMPCFormulations(model);

% Decision variables (nominal)
% u{k} = v_{k-1}, x{k} = z_{k-1}
u = sdpvar(model.nu,secOutputs.N,'full');
x = sdpvar(model.nx,secOutputs.N+1,'full');
% Parameter (measured system state)
X0 = sdpvar(model.nx,1);

% Initialize constraints and objective function
con = [];
obj = 0;

% predefine some variables used in hte MPC definition
AK = model.A + model.B*secOutputs.K;
Xf = secOutputs.params_TVTMPC.Xf;

% ------------------------------ Constraints ------------------------------
% Eq (3.11b): x - z0 \in X_f  =>  r_i^T (x - z0) <= 1   (normalized)
% Here: x is measured state X0, z0 is x(:,1)
con = [con, Xf.A*(X0 - x(:,1)) <= Xf.b];

% Optional: keep measured state inside X bounds (not required by paper)
con = [con, model.x.min <= X0 <= model.x.max];

% Optional: keep nominal initial state z0 inside X (not required by paper)
con = [con, model.x.min <= x(:,1) <= model.x.max];

% Stage loop: paper k = 0..N-1 corresponds to MATLAB idx k = 1...N with kp=k-1
for k = 1:secOutputs.N
    % Eq(3.11a): nominal dynamics
    con = [con, x(:,k+1) == model.A*x(:,k) + model.B*u(:,k)];

    % ---------- Stage constraint tightening ----------
    % For paper index kp = k-1, cross-section is 
    % S_kp = MINKOWSKI_SUM_{j=0}^{kp-1} AK^j W
    % => for k=1 (kp=0): NO subtraction
    % => for k=2 (kp=1): subtract W
    % => for k=3 (kp=2): subtract W + AK W, etc.
    Xset = secOutputs.Px;
    for j = 0:(k-2)
        Xset = Xset - (AK^j)*secOutputs.Pw;
    end

    % Constraint: z_kp + AK^kp (X0 - z0) \in X - S_kp
    % Here: z_kp is x{k}, z0 is x{1}, kp = k-1
    con = [con, Xset.A*(x(:,k) + (AK^(k-1))*(X0 - x(:,1))) <= Xset.b];

    % Input tightening: u = v + K(.) + K*s~, with s~ in S_kp
    Uset = secOutputs.Pu;
    for j = 0:(k-2)
        Uset = Uset - secOutputs.K*(AK^j)*secOutputs.Pw;
    end
    con = [con, Uset.A*(u(:,k) + secOutputs.K*(AK^(k-1))*(X0 - x(:,1))) ...
        <= Uset.b];
end

% % Delta-U constraints
% if MPCforms.deltaucon
%     uprev = sdpvar(model.nu,1,'full');
%     Pdu = secOutputs.Pdu;
%     for k = 1:secOutputs.N
%         if k == 1
%             con = [con, Pdu.A*(u(:,k) - uprev) <= Pdu.b];
%         else
%             con = [con, Pdu.A*(u(:,k) - u(:,k-1)) <= Pdu.b];
%         end
%     end
% end
% Delta-U constraints
if MPCforms.deltaucon
    uprev = sdpvar(model.nu,1,'full');
    Pdu = secOutputs.Pdu;
    for k = 1:secOutputs.N
        if k == 1
            Utube = u(:,k) + secOutputs.K*( X0 - x(:,k) );
            con = [con, Pdu.A*(Utube - uprev) <= Pdu.b];
        else
%             Utube = u(:,k) + secOutputs.K*( x(:,k-1) - x(:,k) );
%             con = [con, Pdu.A*(Utube - Utubeprev) <= Pdu.b];
        end
%         Utubeprev = Utube;
    end
end
% -------------------------------------------------------------------------



% --------------------------- Objective function --------------------------
for k = 1:secOutputs.N
    % Cost: sum z_k^T Q z_k + v_k^T R v_k
    obj = obj + x(:,k)'*model.x.penalty.weight*x(:,k) + ...
                u(:,k)'*model.u.penalty.weight*u(:,k);
end
% ------------------------------------------------------------------------



% --------------------------- Terminal obj-cons ---------------------------
if options.LQRstability
    model = ComputeLQRSet(model);
end

% Terminal constraint
if isprop(model.x,'terminalSet')
    % Eq(3.11d): z_N + AK^N (X0 - z0) \in Xf - S_N
    XN = secOutputs.params_TVTMPC.Xf;
    for k = 0:(secOutputs.N-2)
        XN = XN - (AK^k)*secOutputs.Pw;
    end
    con = [con, XN.A*( x(:,end) + (AK^secOutputs.N)*(X0 - x(:,1)) ) <= XN.b];
end

% Terminal penalty
if isprop(model.x,'terminalPenalty')
    obj = obj + x(:,end)'*model.x.terminalPenalty.weight*x(:,end);
end
% ------------------------------------------------------------------------



% ---------------------------- MPC definition -----------------------------
% Optimization settings - silent verbose mode
opt_settings = sdpsettings('verbose',0);

p_required = X0;
paramtext = 'parameters = [X0';
if any([MPCforms.deltaucon,MPCforms.deltaupen])
    p_required = [p_required; uprev];
    paramtext = ';uprev';
end
paramtext = [paramtext, ']'];

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
OptimProblem.paramtext = paramtext;
% Save all results
secOutputs.OptimProblem = OptimProblem;
end



function  [form] = DeterminTMPCFormulations(model)
% Available MPC forms
form.deltaucon = 0;
form.deltaupen = 0;
form.ytrack = 0;
form.xtrack = 0;
form.xterminalPenalty = 0;
form.ublock = 0;
form.softcon = 0;

if  any([isprop(model.u,'deltaMin'),isprop(model.u,'deltaMax')])
    form.deltaucon = 1;
end

if isprop(model.u,'deltaPenalty')
    form.deltaupen = 1;
end

if isprop(model.x,'reference')
    form.xtrack = 1;
end

if isprop(model.y,'reference')
    form.ytrack = 1;
end

if isprop(model.x,'terminalPenalty')
    form.xterminalPenalty = 1;
end

if isprop(model.u,'block')
    form.ublock = 1;
end

con1 = [isprop(model.x,'softMin'),isprop(model.x,'softMax')];
con2 = [isprop(model.u,'softMin'),isprop(model.u,'softMax')];
con3 = [isprop(model.y,'softMin'),isprop(model.y,'softMax')];
if  any([con1; con2; con3])
    form.softcon = 1;
end
end


function [] = CheckFeasibleSetup(model,options)

form = DeterminTMPCFormulations(model);

switch options.TubeType
    case {'explicit','implicit','timevarying'}
        % Errors
        if form.deltaucon && (options.solType == 0)
            error('MPTplus does not support Tube MPC with delta-u constraints with solType = 0 (use solType = 1), yet.')
        elseif form.deltaupen
            error('MPTplus does not support Tube MPC with delta-u penalization, yet.')
        elseif form.ytrack
            error('MPTplus does not support Tube MPC with y-tracking, yet.')
        elseif form.xtrack
            error('MPTplus does not support Tube MPC with x-tracking, yet.')
        elseif form.ublock
            error('MPTplus does not support Tube MPC with move-blocking, yet.')
        elseif form.softcon
            error('MPTplus does not support Tube MPC with soft constraints, yet.')
        end

        % Warnings
        if form.xterminalPenalty
            warning('MPTplus: Terminal penalty was selected by the user. Check if it was desinged for the nominal system (not the uncertain).')
        end
        if any([options.Nz, all(~isnan(options.Kz),'all'), all(~isnan(options.Pz),'all')])
            warning('MPTplus: Multiple terminal steps "Nz" / controller "Kz" / or penalty "Pz" are not supported in this tube-type. (They do not have effect on the MPC.)')
        end
        if form.deltaucon
             warning('MPTplus: Note that the delta-u constraints are enforced only for the nominal control actions u_opt (not the RHC: Utube = u_opt + K*( x - x_opt )).')
        end
    case 'elastic'
                % Errors
        if form.deltaucon && (options.solType == 0)
            error('MPTplus does not support Tube MPC with delta-u constraints with solType = 0 (use solType = 1), yet.')
        elseif form.deltaupen
            error('MPTplus does not support Tube MPC with delta-u penalization, yet.')
        elseif form.ytrack
            error('MPTplus does not support Tube MPC with y-tracking, yet.')
        elseif form.xtrack
            error('MPTplus does not support Tube MPC with x-tracking, yet.')
        elseif form.ublock
            error('MPTplus does not support Tube MPC with move-blocking, yet.')
        elseif form.softcon
            error('MPTplus does not support Tube MPC with soft constraints, yet.')
        end

        % Warnings
        if form.xterminalPenalty
            warning('MPTplus: Terminal penalty was selected by the user. Check if it was desinged for the nominal system (not the uncertain).')
        end
end
end








