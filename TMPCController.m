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
%      Tube_MaxIter       - maximal number of itterations (default = 1e3)
%      Tube_approx        - {0/1} box approximation of the final Tube
%      solType              - {0/1} defines the retuned optimized values 
%                             (the default value is 0):
%                             0: TMPC returns exact Tube MPC outputs:  
%                                     [u_opt, x_opt] 
%                             1: TMPC returns compact Tube MPC outputs: 
%                                     Utube = u_opt + K*( x - x_opt ),
%                                where [u_opt, x_opt] are the original TMPC
%                                outputs, K is the LQR for the given 
%                                model (given as u = K*x), and x is the  
%                                state measurement.
%     LQRstability          - {0} Does not compute any terminal penalty/set
%                             {1} Computes the terminal penalty/set w.r.t.
%                                 the LTI system (not the ULTI system)
%
% OUTPUTS:
%      TMPC                 - tube MPC policy (MPT3 object)
% ------------------------------------------------------------------------
% % EXAMPLE:
% % LTI system
% model = ULTISystem('A', [1, 1; 0, 1], 'B', [0.5; 1], 'E', [1, 0; 0, 1]);
% model.u.min = [-1]; 
% model.u.max = [ 1];
% model.x.min = [-5; -5]; 
% model.x.max = [ 5;  5];
% model.d.min = [-0.1; -0.1]; 
% model.d.max = [ 0.1;  0.1];
% % Penalty functions
% model.x.penalty = QuadFunction(diag([1, 1]));
% model.u.penalty = QuadFunction(diag([0.01]));
% % Prediction horizon
% N = 2;
% % Include the LQR-based terminal penalty and set for compact control law
% option = {'LQRstability',1, 'solType',1};
% % TMPC controller construction
% TMPC = TMPCController(model,N,option)
% % TMPC evaluation
% x0 = [-5; 2]; % Initial condition
% u = TMPC.evaluate(x0) % MPC evaluation
% ------------------------------------------------------------------------

% REFERENCES:
% [1] https://github.com/holaza/mptplus
% [2] S. Rakovi´c and D. Mayne. A simple tube controller for
% efficient robust model predictive control of constrained
% linear discrete time systems subject to bounded disturbances.
% IFAC Proceedings Volumes, 38(1):241–246,
% 2005. ISSN 1474-6670. doi: https://doi.org/10.3182/20
% 050703-6-CZ-1902.00440. 16th IFAC World Congress.
% [3] Invariant Approximations of the Minimal Robust Positively Invariant 
% Set by S. V. Rakovic, E. C. Kerrigan, K. I. Kouramas, D. Q. Mayne
% IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 50, NO. 3, MARCH 2005
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1406138
% ------------------------------------------------------------------------


%% Checks / initial setup
% Perform checks:
if nargin < 2, error('Ups, not enough inputs were selected'); end
% We do not support parametric uncertainties
if iscell(model.A) || iscell(model.B)
    error('Implemented Tube MPC design method does not support parametric uncertainties (LPV systems), yet.')
end

% set options/default values
ip = inputParser;
% ip.addParameter('eMPC', 0);
ip.addParameter('adjList', 0);
ip.addParameter('Tube_tol', 1e-4);
ip.addParameter('Tube_MaxIter', 1e3);
ip.addParameter('Tube_approx', 0);
ip.addParameter('solType', 1);
ip.addParameter('LQRstability', 0);
if nargin == 3, ip.parse(option{:}); else ip.parse(); end
options = ip.Results;

%% TMPC setup
% Construction of (robust) sets
secOutputs = ConstructSets(model,options);
secOutputs.N = N;

% Define Tube MPC Problem
secOutputs = DefineMPCProblem(secOutputs,model,options);

%% Construction of MPC
% Preconsturct an MPC object
if options.solType
    TMPC = MPCController(model,1);
else
    TMPC = MPCController(model,N);
end
TMPC.construct();

% Construct the Tube (implicit) MPC
TMPC = TubeMPCController( TMPC , secOutputs );

% Overwrite the optimizer based on solType
if options.solType
    TMPC.optimizer = secOutputs.OptimProblem.MPCoptimizers.compact;
else
    TMPC.optimizer = secOutputs.OptimProblem.MPCoptimizers.exact;
end

end


function [secOutputs] = ConstructSets(model,options)
% Function creates sets (constraints) for further processing.
%
% INPUTS:
%           model                   - model of the system (as MPT3 object)
%           options.Tube_tol        - tollerence for the MRPIS
%           options.Tube_MaxIter    - maximal number of itterations for MRPIS
%           options.Tube_approx     - {0,1} request for MRPIS approximation
% OUTPUTS:
%           secOutputs.Tube         - minimal robust positive invariant set
%           secOutputs.Xset         - robustustified state constraint
%           secOutputs.Uset         - robustustified input constraint
%           secOutputs.dUset        - robustustified delta input constraint
%           secOutputs.K            - used LQR policy (as u = K*x)
%           secOutputs.nx           - number of states
%           secOutputs.nu           - number of inputs
%           secOutputs.Xset_full    - state constraint
%           secOutputs.Uset_full    - input constraint

% Bounded additive uncertainty
if isprop(model.d,'setConstraint')
    W = model.d.setConstraint;
else
    W = Polyhedron('lb',model.d.min,'ub',model.d.max);
end
W = model.E*W;
W = Polyhedron('lb',-max(abs(W.V)),'ub',max(abs(W.V)));

% Constraints handling
if isprop(model.x,'setConstraint')
    Xset_full = model.x.setConstraint;
else
    Xset_full = Polyhedron('lb',model.x.min,'ub',model.x.max);
end
if isprop(model.u,'setConstraint')
    Uset_full = model.u.setConstraint;
else
    Uset_full = Polyhedron('lb',model.u.min,'ub',model.u.max);
end
if all([isprop(model.u,'deltaMin'), isprop(model.u,'deltaMax')])
    dUset_full = Polyhedron('lb',model.u.deltaMin ,'ub',model.u.deltaMax );
elseif isprop(model.u,'deltaMin')
    dUset_full = Polyhedron('lb',model.u.deltaMin );
elseif isprop(model.u,'deltaMax')
    dUset_full = Polyhedron('ub',model.u.deltaMax );
else
    dUset_full = {};
end

% Robust positive invariant set
K = model.LQRGain;  % LQR controller as u = K*x
Tube = approx_mRPIS(model.A,model.B,K,W, ...
                options.Tube_tol,options.Tube_MaxIter,options.Tube_approx);
if Tube.isEmptySet
    error('Cannot continue as the minimal positive invariant set is empty!')
end

% State constraints robustification
Xset = Xset_full - Tube;
if Xset.isEmptySet
    error('Cannot continue as the state constraint set is empty!')
end
% Input constraints robustification
UKtube = -K*Tube;
Uset = Uset_full - UKtube;
if Uset.isEmptySet
    error('Cannot continue as the input constraint set is empty!')
end

secOutputs.Tube = Tube;
secOutputs.Xset = Xset;
secOutputs.Uset = Uset;
secOutputs.Xset_full = Xset_full;
secOutputs.Uset_full = Uset_full;
secOutputs.K = K;
secOutputs.nu = model.nu;
secOutputs.nx = model.nx;
end


function [secOutputs] = DefineMPCProblem(secOutputs,model,options)
% Function defines optimization problem of an MPC policy
%
% INPUTS:
%     secOutputs.N          - Prediction horizon
%     secOutputs.Xset_full  - original state consraint
%     secOutputs.Tube       - minimal robust invariant set
%     secOutputs.Xset       - robustustified state constraint
%     secOutputs.Uset       - robustustified input constraint
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
% secOutputs.OptimProblem.MPCoptimizers.compact      - optimizer that retuns 
%                                                      Utube = u_opt + 
%                                                      + K*(x-x_opt) (YALMIP)
% secOutputs.OptimProblem.MPCoptimizers.exact        - optimizer that retuns 
%                                                      [u_opt, x_opt]  (YALMIP)

% Determine MPC formulations that were given by the user
MPCforms = DetermineMPCFormulations(model);


% Define decision variables - YALMIP variables
X0 = sdpvar(model.nx,1);
u = sdpvar(model.nu,secOutputs.N,'full');
x = sdpvar(model.nx,secOutputs.N+1,'full');
% if isprop(model.x,'reference'), xref = sdpvar(model.nx,1,'full'); end
% if  any(isprop(model.u,'deltaPenalty'), ...
%         isprop(model.u,'deltaMin'), ...
%         isprop(model.u,'deltaMax'))
%     uprev = sdpvar(model.nu,1,'full'); 
% end

% Initialize constraints and objective function
con = [];
obj = 0;

% Constraints for OPTIMIZER
con = [con, secOutputs.Xset_full.A*X0 <= secOutputs.Xset_full.b ];

con = [con, secOutputs.Tube.A*( X0 - x(:,1) ) <= secOutputs.Tube.b ];

for k = 1:secOutputs.N
    % Objective - state
    if MPCforms.xtrack
        % TODO: Tube MPC for tracking:
        %  obj = obj + (x(:,k)-xref)'*model.x.penalty.weight*(x(:,k)-xref);
    else
        obj = obj + x(:,k)'*model.x.penalty.weight*x(:,k);
    end
    % Objective - inputs
    if MPCforms.deltaupen
        % TODO: Tube MPC for delta-u formulation:
        % obj = obj + (u(:,k)-uprev)'*model.u.penalty.weight*(u(:,k)-uprev);
    else
        obj = obj + u(:,k)'*model.u.penalty.weight*u(:,k);
    end

    % Nominal LTI model
    con = [con, x(:,k+1) == model.A*x(:,k) + model.B*u(:,k) ];

    % Input and state constraints
    con = [con, secOutputs.Xset.A*x(:,k) <= secOutputs.Xset.b ];
    con = [con, secOutputs.Uset.A*u(:,k) <= secOutputs.Uset.b ];
end

% LQR stability:
% Introduce LQR-based terminal set and terminal penalty evaluated for LTI
% system, not ULTI system.
if options.LQRstability
    model = ComputeLQRSet(model);
end

% Terminal penalty
if isprop(model.x,'terminalPenalty')
    obj = obj + x(:,end)'*model.x.terminalPenalty.weight*x(:,end);
end

% Terminal constraint
if isprop(model.x,'terminalSet')
    con = [con, model.x.terminalSet.A*x(:,end) <= model.x.terminalSet.b ];
end

% Delta-U constraints
if MPCforms.deltaucon
    uprev = sdpvar(model.nu,1,'full');
    dUset = secOutputs.dUset;
    for k = 1:secOutputs.N
        if k == 1
            con = [con, dUset.A*(u(:,k) - uprev) <= dUset.b];
        else
            con = [con, dUset.A*(u(:,k) - u(:,k-1)) <= dUset.b];
        end
    end
end

% Optimization settings - silent verbose mode
opt_settings = sdpsettings('verbose',0);

p_required = X0;

% ------------------------- This is not needed ----------------------------
% OptimProblem.Constraints = con;
% OptimProblem.Objective = obj;
% OptimProblem.Param_Required = p_required;
% OptimProblem.CompactForm_Param_Returned = u(:,1) + secOutputs.K*( X0 - x(:,1));
% % OptimProblem.ExactForm_Param_Returned = [u(:,1); x(:,1)];
% OptimProblem.ExactForm_Param_Returned = [u(:); x(:)];
% OptimProblem.YALMIP_settings = opt_settings;
% -------------------------------------------------------------------------

% Note that the eMPC will have always p_return = [u(:,1); x(:,1)]. If
% required, by solType == 1, then this primal function will be overwritten
% inside of the .toExplicit function
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
    % Import parameters from ULTI system to LTI system
    LTI_model = LTISystem('A', model.A, 'B', model.B);
    LTI_model.x.penalty = QuadFunction( model.x.penalty.weight );
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
    [ Klqr, P ] = dlqr( model.A, model.B, model.x.penalty.weight, model.u.penalty.weight );
    K = -Klqr;
    LTI_model.x.with('terminalPenalty');
    LTI_model.x.terminalPenalty = QuadFunction( P );
    % Export parameters from LTI system to ULTI system
    % Note, the original Termianl Set and Termianl Penalty are overwirtten
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

if k == 1, fprintf(1, strcat(operation,'       :')); end
    
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


function F_alpha = approx_mRPIS(A,B,K,W,tol,MaxIter,Approx)
% Invariant Approximations of the Minimal Robust Positively Invariant Set
% by S. V. Rakovic, E. C. Kerrigan, K. I. Kouramas, D. Q. Mayne
% IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 50, NO. 3, MARCH 2005
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1406138

if nargin < 4, error('Ups, not enough inputs...'); end
if nargin < 5, tol = 1e-4; end
if nargin < 6, MaxIter = 1000; end
if nargin < 7, Approx = 0; end

% Closed-loop/autonomous system matrix 
Acl = A + B*K;

% Problem size
nx = size(Acl, 1); 

% Initial values
s = 0; % "s" is integer going toward zero: s -> 0
alpha = 1e6; % Initialize "alpha"
Ms = 1000;
MSS = zeros(2*nx,MaxIter);

fprintf('Constructing approximated minimal robust invariant set ...\n')
tic,
while( alpha > tol/( tol + Ms ) ) && ( s < MaxIter )

    % Increment "s"
    s = s + 1;

    % Eq(10): (Acl^s)*W \subseteq alpha*W == h_W( (Acl^(s))'*f_i ) <= alpha*g_i
    % Eq(11): alpha(s) = max_{i \in I} ( ( h_W( (Acl^(s))'*f_i ) ) / g_i )
    alpha = max(W.support(Acl^(s)*(W.A)')./W.b);

    % Eq(13): M(s) = max_{j=1,...,n} ( sum_{i=0}^{s-1} ( h_W( (Acl^(i))'*e_j ) ), sum_{i=0}^{s-1} ( h_W( -(Acl^(i))'*e_j ) ) )
    MSS(:,s) = W.support([Acl^(s), -Acl^(s)]);
    mss = sum(MSS(:,1:s),2);
    Ms = max(mss);

    temp = alpha - tol/( tol + Ms );
    if mod(s,10) == 0
        fprintf('\tConvergance to 0: %.5f\n',max(0,temp))
    end
end

% Eq(2): Fs = Minkowski_Sum_{i=0}^{s-1}( Acl^(i)*W ), where F_0 = {0}
Fs = W; 
for i = 1 : (s - 1)
    ShowProgress('Finalizing...         ',i,s-1)
    Fs = Fs + Acl^(i)*W;
    Fs.minVRep();
%     if size(Fs.V,1) > 2e5
%         Fs = Fs.outerApprox;
%     end
end
F_alpha = (1 - alpha)^(-1)*Fs;
compTime = toc;
fprintf('\b (computation time: %.2f, #itterations: %.0f, #halfspaces = %.0f)\n',compTime,s,size(Fs.H,1))

if Approx
    vol = Fs.volume;
    Fs = Fs.outerApprox;
    vol2 = Fs.volume;
    increaseVOl = (vol2-vol)/vol*100;
    fprintf('The Minimal Robust Invariant Set was approximated! (Volume increased by %.0f%%)\n',increaseVOl)
    F_alpha = (1 - alpha)^(-1)*Fs;
end
end 


function  [form] = DetermineMPCFormulations(model)
% Available MPC forms
form.deltaucon = 0;
form.deltaupen = 0;
form.ytrack = 0;
form.xtrack = 0;

if  any([isprop(model.u,'deltaMin'),isprop(model.u,'deltaMax')])
    form.deltaucon = 1;
    fprintf('MPC formulation: Delta-u constraints detected !\n')
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