classdef ExplicitTubeMPCController < EMPCController
    % Class representing explicit tube MPC controllers
    %
    % Constructor:
    %   ctrl = ExplicitTubeMPCController(explicit_controller)
    %

    properties
        TMPCparams
    end

    methods

        function out = getName(obj)
            out = 'Explicit tube MPC controller';
        end

        function obj = ExplicitTubeMPCController( ctrl , secOutputs )
            % Constructor:
            if nargin==0
                return
            elseif ~isa(ctrl, 'EMPCController')
                error('Input must be an EMPCController object.');
            end
            % copy data from the source object
            obj.importUserData(ctrl.copy());

            % Assgin properties of Tube MPC
            if( nargin >= 2 )
                obj.TMPCparams = secOutputs;  % TMPC parameters and solution from YALMIP
            else
                obj.TMPCparams = NaN;
            end
        end


        function [u, feasible, openloop] = evaluate(obj, xinit, varargin)

            [u, feasible, openloop] = obj.evaluate@EMPCController(xinit, varargin{:});

            if feasible
                % Generate OPENLOOP
                if ( length(u) == obj.TMPCparams.nu )
                    %% solType == 1
                    openloop.cost = NaN;
                    openloop.U = NaN;
                    openloop.X = NaN;
                    openloop.Y = NaN;
                    openloop.info = 'OPENLOOP evaluation of explicit Tube MPC is not supported, yet.';
                else
                    %% solType == 0
                    openloop.cost = NaN;
                    openloop.U = NaN;
                    openloop.X = NaN;
                    openloop.Y = NaN;
                    openloop.info = 'OPENLOOP evaluation of explicit Tube MPC is not supported, yet.';
                end
            end
        end


        function [obj] = createAdjList(obj)
            % Creates the adjacency list for the given partition
            % The resulting list can by found in:
            %               obj.optimizer.adj_list

            adj_list = cell(1,obj.nr);
            for k = 1:obj.nr
                ShowProgress('Creating adjacency list ... ', k, obj.nr)
                t = 0;
                temp = cell(1,1);
                for kk = 1:obj.nr
                    if obj.partition.Set(k).isAdjacent(obj.partition.Set(kk))
                        intersection = obj.partition.Set(k) & obj.partition.Set(kk);
                        if size(intersection.V,1) >= obj.partition.Set(k).Dim
                            t = t+1;
                            temp{t,1} = kk;
                        end
                    end
                end
                adj_list{k} = temp;
            end
            % adj_list = MPT_verify_graph(eMPC.partition.Set,adj_list);
            obj.optimizer.setInternal('adj_list', adj_list);
        end

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
                % Parse inputs
                mpc = obj;
                % Extract parameters
                nx = mpc.TMPCparams.nx;
                nu = mpc.TMPCparams.nu;
                K = mpc.TMPCparams.K;
                W = mpc.TMPCparams.Wset;
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
                    % Uncertain (noisy) LTI system: X(k+1) = A*X(k) + B*U(k)
                    Xdata(:,k+1) = A*X + B*( Utube ) + w(:,k);
                    % Output data
                    Unominal = [Unominal, U_opt];
                    Xnominal = [Xnominal, X_opt];
                    Udata = [Udata, Utube];
                    cost = cost + X'*Q*X + Utube'*R*Utube; % Quadratic cost
                end
                % Assign outputs
                ClosedLoopData.U = Udata;
                ClosedLoopData.X = Xdata;
                ClosedLoopData.cost = cost;
                ClosedLoopData.Unominal = Unominal;
                ClosedLoopData.Xnominal = Xnominal;
                ClosedLoopData.K = K;
            else
                [ ClosedLoopData ] = obj.simulate@EMPCController(xinit, Nsim, varargin{:});
            end % if ( solType == 0 )
        end % function

        function [] = indexPartition( obj, indices2depict, figureHandle )
            % Function indexPartition plot parition and  depicts indices of each
            % critical region.
            obj.partition.plot
            depictIndices( obj, [1 : obj.nr], gcf )
        end

    end
end
