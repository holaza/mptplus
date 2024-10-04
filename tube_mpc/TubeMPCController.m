classdef TubeMPCController < MPCController
    % Class representing on-line tube MPC controllers
    %
    % Constructor:
    %   ctrl = TMPCController(model, N)
    %          model: prediction model
    %              N: prediction horizon

    properties
        TMPCparams
    end

    methods

        function out = getName(obj)
            out = 'Tube MPC controller';
        end

        function obj = TubeMPCController( ctrl, secOutputs )
            % Constructor:
            %   ctrl = Controller(model, N)
            %         model: prediction model
            %             N: prediction horizon

            if nargin==0
                return
            elseif ~isa(ctrl, 'MPCController')
                error('Input must be an MPCController object.');
            end

            obj.importUserData(ctrl.copy());

            % Assgin properties of Tube MPC
            if( nargin >= 2 )
                obj.TMPCparams = secOutputs;  % TMPC parameters and solution from YALMIP
            else
                obj.TMPCparams = NaN;
            end
        end

        function EMPC = toExplicit(obj)
            % Extract all required parameters from the obj
            option = obj.TMPCparams.OptimProblem.options;
            model = obj.model;
            K = obj.TMPCparams.K;

            % Compute the explicit controller
            mpsol = obj.TMPCparams.OptimProblem.MPCproblem.solve;
            Tempc = mpsol.xopt;
            if mpsol.exitflag ~= 1,error('Ups, the problem is infeasible!');end

            % A) Convert the feedback law to the implicit form
            if option.solType
                % A) Update the feedback policy
                % Primal function: u =  Fi*x + gi
                % Tube MPC primal function: Utube = u_opt + K*( x - x_opt ), where:
                %   x_opt - computed initial state of the prediction (encoded in u)
                %   u_opt - computed control inputs
                %   K     - assumed LQR policy (from the TMPC setup)
                %   u     - u = [u_opt; x_opt]
                % ------------------------------------------------------
                % Primal function: Utube = (Fu_opt + K - K*Fx_opt)*x + (gu_opt-K*gx_opt)
                indx_x_opt = model.nu+1:model.nu+model.nx;
                indx_u_opt = 1:model.nu;
                % Extending the LQR gain to cover also all other parameters
                KK = [K, zeros(1,Tempc.Set(1).Dim - model.nx)];
                for k = 1:Tempc.Num
                    F = Tempc.Set(k).getFunction('primal').F;
                    g = Tempc.Set(k).getFunction('primal').g;
                    % F = (Fu_opt + K - K*Fx_opt)
                    Fnew = F(indx_u_opt,:) + KK - K*F(indx_x_opt,:);
                    % g = (gu_opt - K*gx_opt)
                    gnew = g(indx_u_opt,:) - K*g(indx_x_opt,:);
                    % Overwrtite the old 'primal'function
                    Tempc.Set(k).addFunction(AffFunction(Fnew,gnew), 'primal');
                end
            end

            % B) Convert the solutino into MPT3 object
            EMPC = EMPCController(Tempc);

            % C) Recover lost information
            EMPC.N = 1; % To satisfy the sanity check of function .EVALUATE (N*nu)
            EMPC.model = model;

            % D)
            EMPC = ExplicitTubeMPCController( EMPC , obj.TMPCparams );

            % E) Change the number of inputs (needed for '.evaluate')
            if ~option.solType
                EMPC.nu = model.nu + model.nx;
            end

            % F) Add list of xinitFormat (for .evaluate())
            EMPC.xinitFormat = obj.xinitFormat;
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

            error(validate_vector(xinit, obj.nx, 'initial state'));

            % assemble the vector of initial conditions. Include any
            % variables that were declared by filters as those which need
            % proper initialization.
            xinit = obj.parse_xinit(xinit, varargin{:});

            % use a pre-constructed optimizer object for faster
            % evaluation
            [U, status] = obj.optimizer{xinit};
            % these statuses indicate a feasible solution
            % (from "help yalmiperror"):
            %  0: Successfully solved
            %  3: Maximum #iterations or time-limit exceeded
            %  4: Numerical problems
            %  5: Lack of progress
            feasible = ismember(status, [0, 3, 4, 5]);
            if ~feasible
                J = NaN;
                u = NaN(obj.nu, 1);
                U = NaN(size(U));
            elseif isempty(U)
                error('MPTplus: Ooops, something went wrong. Please report to the contact e-mail address.');
            else
                J = U(1);

                solType = obj.TMPCparams.OptimProblem.options.solType;
                if ( solType == 1 )
                    % Return compact control law: u(t) = u0 + K*(x(t) - x0)
                    u = U(2 : obj.TMPCparams.nu+1);
                elseif ( solType == 0 )
                    switch obj.TMPCparams.OptimProblem.options.TubeType
                        case {'explicit','implicit'}
                            % Return exact control law: [u0, x0]

                            % U_exact = reshape( U(2 : nu*N + 1 ), [nu, N] )
                            U_exact = reshape ( U( 2 : obj.TMPCparams.nu*obj.TMPCparams.N + 1 ), [ obj.TMPCparams.nu, obj.TMPCparams.N ] );

                            % X_exact = reshape( U(nu*N + 2 : end), [nx, N+1] )
                            X_exact = reshape ( U( obj.TMPCparams.nu*obj.TMPCparams.N + 2 : end) , [ obj.TMPCparams.nx, obj.TMPCparams.N+1 ] );
                            u = [ U_exact( : , 1 ); X_exact( : , 1 ) ]; % return 1st contol step: [u0, x0]
                        case 'elastic'
                            % Return exact control law: [u0, x0, sigma0]

                            % U_exact = reshape( U(2 : nu*N + 1 ), [nu, N] )
                            start = 2;
                            fin = start + obj.TMPCparams.nu*obj.TMPCparams.N - 1;
                            U_exact = reshape ( U( start : fin ), [ obj.TMPCparams.nu, obj.TMPCparams.N ] );

                            % X_exact = reshape( U(nu*N + 2 : nu*N + 2 + nx*(N+1) - 1), [nx, N+1] )
                            start = fin + 1;
                            fin = start + obj.TMPCparams.nx*(obj.TMPCparams.N + 1) -1;
                            X_exact = reshape ( U( start : fin) , [ obj.TMPCparams.nx, obj.TMPCparams.N+1 ] );
                            
                            % sigma_exact = reshape( U(nu*N + 2 + nx*(N+1) : nu*N + 2 + nx*(N+1) + qs*N  - 1), [qs, N] )
                            start = fin + 1;
                            fin = start + obj.TMPCparams.elasticparams.qs*(obj.TMPCparams.N) - 1;
                            sigma_exact = reshape ( U( start : fin) , [ obj.TMPCparams.elasticparams.qs, obj.TMPCparams.N ] );
                            
                            u = [ U_exact( : , 1 ); X_exact( : , 1 ); sigma_exact(:,1) ]; % return 1st contol step: [u0, x0, sigma0]
                        otherwise
                            error('MPTplus: Incorrect TubeType found during control action evaluation!')
                    end
                else
                    error(sprintf('MPTplus: Unexpected value of solType: %d',solType))
                end
            end

            if ( nargout == 3 )
                solType = obj.TMPCparams.OptimProblem.options.solType;
                if ( solType == 1 )
                    openloop.cost = NaN;
                    openloop.U = NaN;
                    openloop.X = NaN;
                    openloop.Y = NaN;
                    openloop.info = 'OPENLOOP evaluation of Tube MPC is under cosntruction.';
                elseif ( solType == 0 )
                    switch obj.TMPCparams.OptimProblem.options.TubeType
                        case {'explicit','implicit'}
                            % U' = [ cost, [u(0)',...,u(N-1)'], [x(0)',...,x(N)'] ]
                            openloop.cost = J;
                            openloop.U = U_exact;
                            openloop.X = X_exact;
                            openloop.Y = NaN;
                            openloop.info = 'OPENLOOP evaluation of Tube MPC is under cosntruction.';
                        case 'elastic'
                            % U' = [ cost, [u(0)',...,u(N-1)'], [x(0)',...,x(N)'], [sigma(0),...,sigma(N-1)] ]
                            openloop.cost = J;
                            openloop.U = U_exact;
                            openloop.X = X_exact;
                            openloop.sigma = sigma_exact;
                            openloop.Y = NaN;
                            openloop.info = 'OPENLOOP evaluation of Tube MPC is under cosntruction.';
                        otherwise
                            error('MPTplus: Incorrect TubeType found during control action evaluation!')
                    end
                else
                    error(sprintf('MPTplus: Unexpected value of solType: %d',solType))
                end
            end
        end % function

        function [ ClosedLoopData ] = simulate(obj, xinit, Nsim, varargin)
            % [ ClosedLoopData ] = simulate(obj, xinit, Nsim, varargin)
            % If Tube MPC design handles the expanded control law (solType==0),
            % then the function .SIMULATE evaluates the corresponding
            % simulation of the closed-loop control.
            % Otherwise, a conventional .SIMULATE of MPT toolbox is called.
            if( nargin < 3 ) % Check number of inputs
                error(' The closed-loop simulation requires at least 3 inputs: MPC controller, x0, Nsim!')
            end
            % Check type (solType) of Tube MPC controller: "u_tube" vs "[u0, x0]"
            solType = obj.TMPCparams.OptimProblem.options.solType;
            if ( solType == 0 )
                tubeType = obj.TMPCparams.OptimProblem.options.TubeType;
                % Parse inputs
                mpc = obj;
                % Extract parameters
                nx = mpc.TMPCparams.nx;
                nu = mpc.TMPCparams.nu;
                K = mpc.TMPCparams.K;
                W = mpc.TMPCparams.Pw;
                A = mpc.model.A;
                B = mpc.model.B;
                Q = mpc.model.x.penalty.H;
                R = mpc.model.u.penalty.H;
                % Sequence of additive disturbances
                if( isempty(varargin) == 1 )
                    % Generate random sequence of disturbance
                    for k = 1 : Nsim
                        w(:,k) = W.randomPoint;
                    end
                else
                    % Load sequence of disturbance from function inputs
                    w = varargin{1};
                end
                % Initialization of the closed-loop control outputs
                X = xinit; % initial system state
                Unominal = []; % initial nominal control action
                Xnominal = [ xinit ]; % initial nominal system states
                Udata = []; % initial closed-loop control action
                Xdata = [ X ]; % initial closed-loop system states
                if strcmpi(tubeType,'elastic'), SIGMAnominal = []; end
                cost = 0;
                for k = 1 : Nsim
                    % Current system states
                    X = Xdata(:,k);
                    % Current control action
                    [ UX_opt ] = mpc.evaluate(X);
                    % Extract variables
                    U_opt = UX_opt(1:nu);
                    X_opt = UX_opt(nu+1 : nu+nx);
                    Utube = U_opt + K*( X - X_opt );
                    if strcmpi(tubeType,'elastic') 
                        Sigma = UX_opt(nu+nx+1:nu+nx + mpc.TMPCparams.elasticparams.qs); 
                        Utube = U_opt + mpc.TMPCparams.elasticparams.Ka(Sigma)*( X - X_opt );
                    end
                    % Uncertain (noisy) LTI system: X(k+1) = A*X(k) + B*U(k)
                    Xdata(:,k+1) = A*X + B*( Utube ) + w(:,k);
                    % Output data
                    Unominal = [Unominal, U_opt];
                    Xnominal = [Xnominal, X_opt];
                    Udata = [Udata, Utube];
                    if strcmpi(tubeType,'elastic')
                        SIGMAnominal = [SIGMAnominal,Sigma]; 
                    end
                    cost = cost + X'*Q*X + Utube'*R*Utube; % Quadratic cost
                end
                % Assign outputs
                ClosedLoopData.U = Udata;
                ClosedLoopData.X = Xdata;
                ClosedLoopData.W = w;
                ClosedLoopData.cost = cost;
                ClosedLoopData.Unominal = Unominal;
                ClosedLoopData.Xnominal = Xnominal;
                ClosedLoopData.K = K;
                if strcmpi(tubeType,'elastic') 
                    ClosedLoopData.K = mpc.TMPCparams.elasticparams.Ka();
                    ClosedLoopData.SIGMAnominal = SIGMAnominal; 
                end
            else
                %% if ( solType == 1 ) 
                [ ClosedLoopData ] = obj.simulate@MPCController(xinit, Nsim, varargin{:});
                
                % Overwrite original MPT3 "cost" (NaN) with the closed loop cost 
                cost = 0;
                Q = obj.model.x.penalty.H;
                R = obj.model.u.penalty.H;
                for k = 1 : Nsim
                    % Current system states
                    X = ClosedLoopData.X(:,k);
                    % Current control action
                    U = ClosedLoopData.U(:,k);
                    cost = cost + X'*Q*X + U'*R*U; % Quadratic cost
                end
                ClosedLoopData.cost = cost; 

            end % if ( solType == 0 )
        end % function
    end

end
