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
            % TODO: Replace by a stable construction of the Explicit Tube MPC
            error('MPTplus: Construction of the the Explicit Tube MPC is not supported, yet.')
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
                J = U(1); % cost

                solType = obj.TMPCparams.OptimProblem.options.solType;
                if ( solType == 1 )
                    % Return compact control law: u(t) = u0 + K*(x(t) - x0)
                    u = U(2 : obj.TMPCparams.nu+1);
                elseif ( solType == 0 )
                    % Return exact control law: [u0, x0]

                    % U_exact = reshape( U(2 : nu*N + 1 ), [nu, N] )
                    U_exact = reshape ( U( 2 : obj.TMPCparams.nu*obj.TMPCparams.N + 1 ), [ obj.TMPCparams.nu, obj.TMPCparams.N ] );

                    % X_exact = reshape( U(nu*N + 2 : end), [nx, N+1] )
                    X_exact = reshape ( U( obj.TMPCparams.nu*obj.TMPCparams.N + 2 : end) , [ obj.TMPCparams.nx, obj.TMPCparams.N+1 ] );
                    u = [ U_exact( : , 1 ); X_exact( : , 1 ) ]; % return 1st contol step: [u0, x0]
                else
                    error(sprintf('MPTplus: Unexpected value of solType: %d',solType))
                end
			end

			if nargout==3
                if ( solType == 1 )
                    openloop.cost = NaN;
                    openloop.U = NaN;
                    openloop.X = NaN;
                    openloop.Y = NaN;
                    openloop.info = 'OPENLOOP evaluation of Tube MPC is under cosntruction.';
                elseif ( solType == 0 )
                    % U' = [ cost, [u(0)',...,u(N-1)'], [x(0)',...,x(N)'] ]
                    openloop.cost = J;
                    openloop.U = reshape( U( 2:( obj.TMPCparams.N * obj.TMPCparams.nu )+1) , [ obj.TMPCparams.nu, obj.TMPCparams.N ] );
                    openloop.X = reshape( U( ( obj.TMPCparams.N * obj.TMPCparams.nu )+2 : end ) , [ obj.TMPCparams.nx, obj.TMPCparams.N+1 ] );
                    openloop.Y = NaN;
                    openloop.info = 'OPENLOOP evaluation of Tube MPC is under cosntruction.'
                else
                    error(sprintf('MPTplus: Unexpected value of solType: %d',solType))
                end
			end
        end
    end
end
