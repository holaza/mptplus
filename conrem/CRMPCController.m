classdef CRMPCController < MPCController



    properties
        CRparams        % stored parameters of the controller
        InActiveCons    % denotes which con is still inactive
                        % (i.e. cons(1:InActiveCons) are inactive)
        Jprev           % cost value of the previous itteration
        ConStat         % denotes if the controller is fully initialized:
                        %     (0) has removed cons
                        %     (1) is fully initialized (all cons)
       
    end

    properties (Access = private)
        requestedFormat = [];
        estimator = [];
    end

    methods

        function out = getName(obj)
            out = 'Contraints removal MPC controller';
        end

        function obj = CRMPCController( ctrl, params )
            % Constructor:

            if nargin==0
                return
            elseif ~isa(ctrl, 'MPCController')
                error('Input must be an MPCController object.');
            end

            % copy essential information from the ctrl
            obj.importUserData(ctrl.copy());
            obj.requestedFormat = params.optprob.requestedFormat;

            % Assgin properties of Constraints Removal MPC
            if( nargin >= 2 )
                obj.CRparams = params;  % parameters and solution from YALMIP
                obj.InitializeCRMPC(obj.CRparams.options.initializeType)
            else
                error('Not enough inputs ...')
            end
        end

        function InitializeCRMPC(obj, option)
            if nargin < 2, option = 0; end
            if option
                % initialization with removed "infeasible" constraints
                temp = sum(obj.CRparams.sigma == inf);
                obj.InActiveCons = obj.CRparams.optprob.nineq - temp;
                fprintf('Controller was initialized, %.0f/%.0f inequality constraints were removed due to redundancy.\n',temp,obj.CRparams.optprob.nineq)
            else
                obj.InActiveCons = obj.CRparams.optprob.nineq;
                fprintf('Controller was initialized with %.0f/%.0f inequality constraints.\n',obj.CRparams.optprob.nineq,obj.CRparams.optprob.nineq)
            end
            obj.Jprev = inf;
            obj.ConStat = 1;
        end

        function [answr] = isFullyInitialized(obj)
            answr = obj.ConStat == 1;
        end

        function [u, feasible, openloop] = evaluate(obj, xinit, varargin)
            % Solves the MPC optimization problem for a given initial
            % condition
            %
            % u = controller.evaluate(x0) solves the MPC optimization
            % problem using x0 as the initial condition and returns the
            % first element of the optimal sequence of control inputs. If
            % the problem is infeasible, "u" will be NaN.
            %
            % [u, feasible] = controller.evaluate(x0) also returns a
            % boolean flag indicating whether the optimization problem had
            % a feasible solution.
            %
            % [u, feasible, openloop] = controller.evaluate(x0) also
            % returns the open-loop predictions of states, inputs and
            % outputs in "openloop.X", "openloop.Y", and "openloop.Y",
            % respectively. Value of the optimal cost is returned in
            % openloop.cost.
            %
            % u = controller.evaluate(x0, 'x.reference', ref, 'u.prev', u0)
            % includes "ref" and "u0" into the vector of initial
            % conditions.

            % make sure we have the prediction horizon available
            error(obj.assert_controllerparams_defined);
          
			% assemble the vector of initial conditions. Include any
			% variables that were declared by filters as those which need
			% proper initialization.
            if ~isequal(size(xinit), [obj.xinitFormat.n_xinit 1])
                xinit = obj.parse_xinit(xinit, varargin{:});
            end

            % -------------------- Constraints removal --------------------
            % Remove constraints (if CR is NOT set to fixed)
            if ~obj.CRparams.options.CRfixed
                tic;
                % find constraints to remove
                temp = sum(obj.CRparams.sigma(1:obj.InActiveCons) > obj.Jprev);
                obj.InActiveCons = obj.InActiveCons - temp;
                % report that the contorller is no longer initialized (i.e. is modified)
                if obj.CRparams.optprob.nineq > obj.InActiveCons
                    obj.ConStat = 0;
                end
                overheadtime = toc;
            else
                overheadtime = 0;
            end
            % -------------------- Constraints removal --------------------


            % Solve the problem
            [U, info] = timed_quadprog(obj,xinit);

            % Quadporg status:
            %     All algorithms:
            %       1  First order optimality conditions satisfied.
            %       0  Maximum number of iterations exceeded.
            %      -2  No feasible point found.
            %      -3  Problem is unbounded.
            %     Interior-point-convex:
            %       2  Solver stalled at feasible point.
            %      -6  Non-convex problem detected.
            %      -8  Unable to compute step direction; no further progress can be made.
            %     Trust-region-reflective:
            %       3  Change in objective function too small.
            %      -4  Current search direction is not a descent direction; no further
            %          progress can be made.
            %     Active-set:
            %      -6  Non-convex problem detected.
            feasible = ismember(info.problem, [0, 1, 2, 3]);

            if ~feasible
                u = NaN(obj.nu, 1);
            elseif isempty(U)
                error('MPTplus: Ooops, something went wrong. Please report to the contact e-mail address.');
            else
                u = U(obj.requestedFormat.indices.u(1:obj.nu));
                obj.Jprev = U(1);
            end

            if ( nargout == 3 )
                indexu = strcmp(obj.requestedFormat.names, 'u');
                indexx = strcmp(obj.requestedFormat.names, 'x');
                if ~feasible
                    UU = NaN(obj.nu*obj.N,1);
                    X = NaN(obj.nx*(obj.N+1),1);
                    J = NaN;
                    Z = NaN(obj.CRparams.optprob.nz,1);
                else
                    % Need to reconstruct these outputs: u, feasible, openloop
                    UU = U(obj.requestedFormat.indices.u);
                    X = U(obj.requestedFormat.indices.x);
                    J = U(obj.requestedFormat.indices.J);
                    Z = U(2:end);
                end
                openloop.cost = J;
                openloop.U = reshape(UU,obj.requestedFormat.dims{indexu});
                openloop.X = reshape(X,obj.requestedFormat.dims{indexx});
                openloop.Z = Z;

                % Append Y if needed
                if any(strcmp(obj.requestedFormat.names, 'y'))
                    indexy = strcmp(obj.requestedFormat.names, 'y');
                    if ~feasible
                        Y = NaN(obj.model.ny*obj.N,1);
                    else
                        Y = U(obj.requestedFormat.indices.y);
                    end
                    openloop.Y = reshape(Y,obj.requestedFormat.dims{indexy});
                end
                
                % add references if needed
                con1 = obj.CRparams.optprob.MPCFormulation.xtracking;
                con2 = obj.CRparams.optprob.MPCFormulation.ytracking;
                if con1
                    references.x = U(obj.CRparams.optprob.requestedFormat.zindices.xref,:);
                end
                if con2
                    references.y = U(obj.CRparams.optprob.requestedFormat.zindices.yref,:);
                end
                if con1 || con2
                    openloop.references = references;
                end

                % Append additional information
                openloop.solvertime = info.solvertime;
                openloop.ncons = obj.InActiveCons;
                openloop.overheadtime = overheadtime;
            end
        end



function out = simulate(obj, x0, N_sim, varargin)
    % Simulate the closed-loop system using the prediction model
    %
    %   data = controller.simulate(x0, N_sim)
    %
    % Inputs:
    %   controller: a controller object
    %   x0: initial point for the simulation
    %   N_sim: number of simulation steps
    %
    % Outputs:
    %   data: structure containing closed-loop profiles of states,
    %         inputs, outputs, and the cost function
    %
    % See ClosedLoop/simulate for more information.

    out = ClosedLoopCR(obj, obj.model).simulate(x0, N_sim, varargin{:});
end

    end


    methods (Access = protected)
        function [U, info] = timed_quadprog(obj,x0)
            % The main pourpose of this function is to measure evaluation
            % time while executing quadprog only once.
            % TODO: Find a better way to measure evaluatikon time!

            % Store the result here
            persistent mptplus_zopt
            persistent mptplus_fval
            persistent mptplus_exitflag
            persistent mptplus_solvtime
            t_helper = [];
            
            if exist('mosekopt.m','file')
                options = [];
            else
                options = optimoptions('quadprog', 'Display', 'none');
            end

            % Wrapper function for timeit
            function quadprog_only()
                [mptplus_zopt, mptplus_fval, mptplus_exitflag,mptplus_solvtime] = quadprog(2*obj.CRparams.optprob.Q, ...
                    obj.CRparams.optprob.b, ...
                    obj.CRparams.optprob.Aineq(1:obj.InActiveCons,:), ...
                    obj.CRparams.optprob.bineq(1:obj.InActiveCons), ...
                    [obj.CRparams.optprob.Aeq; obj.CRparams.optprob.Ainit], ...
                    [obj.CRparams.optprob.beq; x0], ...
                    [], [], [], options);
            end

            % Measure time
            switch obj.CRparams.options.SolveType
                case 0
                    t_helper = clock;
                    quadprog_only();
                    execution_time = etime(clock,t_helper);
                    % try to get the build-in quadprog timed eval time
                    % (available from Matlab ~2023)
                    try
                        execution_time = mptplus_solvtime.time;
                    end 
                case 1
                    execution_time = timeit(@quadprog_only);
                otherwise
                    error('MPT+: Ups, something went wrong!')
            end
            

            % Retrieve the result
            U = [mptplus_fval; mptplus_zopt];
            info.problem = mptplus_exitflag;
            info.solvertime = execution_time;
            info.info = mptplus_solvtime;
        end
    end
end
