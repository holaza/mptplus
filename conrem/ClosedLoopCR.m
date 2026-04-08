classdef ClosedLoopCR < ClosedLoop

    methods
        function obj = ClosedLoopCR(system, controller)
            obj@ClosedLoop(system, controller);
        end

        function out = simulate(obj, x0, N_sim, varargin)
            %
            % Simulates the closed-loop system for a given number of steps
            %

            narginchk(3, Inf);
            error(validate_vector(x0, obj.system.nx, 'initial state'));

            % internal helper to derimine if "string" ends up with "part"
            function out = endswith(string, part)
                out = length(string)>=length(part) && ...
                    isequal(string(end-length(part)+1:end), part);
            end

            % determine whether we have free references. if so, allow the
            % user to specify their respective full closed-loop profiles
            references = struct('value', {}, 'position', {});
            if nargin>3
                if mod(length(varargin), 2)~=0
                    error('Options must come in key/value pairs.');
                end
                % find position of reference options in varargin
                ref_positions = find(cellfun(@(z) endswith(z, '.reference'), varargin));
                for i = 1:length(ref_positions)
                    % references found, store their values and position of
                    % the value in varargin such that we can updated them
                    % later
                    ref_name = varargin{ref_positions(i)};
                    ref_value = varargin{ref_positions(i)+1};
                    % validate dimensions: number of columns of ref_value
                    % must either be 1 or Nsim
                    if size(ref_value, 2)==1
                        % expand throughout N_sim
                        ref_value = repmat(ref_value, 1, N_sim);
                    elseif size(ref_value, 2)~=N_sim
                        error('"%s" must have either 1 or %d columns.', ...
                            ref_name, N_sim);
                    end
                    % store the reference
                    ref.value = ref_value;
                    % position of the reference value in varargin
                    ref.position = ref_positions(i)+1;
                    references(end+1) = ref;
                end
            end

            include_disturbance = isa(obj.system, 'ULTISystem');

            X = x0(:); U = []; Y = []; D = []; J = [];
            % -------------------------------------------------------------
            pred_ctrl = []; % open-loop controller predictions at each step
            SolverTime = zeros(N_sim,1);                                    % updated
            ncons = zeros(N_sim,1);                                         % updated
            Zdata = nan(obj.controller.CRparams.optprob.variables.dim.nz,N_sim);        % updated
            % -------------------------------------------------------------
            for k = 1:N_sim

                if k>1
                    % update u.previous and y.previous
                    u_position = find(cellfun(@(z) isequal(z, 'u.previous'), varargin));
                    if ~isempty(u_position)
                        varargin{u_position+1} = u;
                    end
                    y_position = find(cellfun(@(z) isequal(z, 'y.previous'), varargin));
                    if ~isempty(y_position)
                        varargin{y_position+1} = y;
                    end
                    % TODO: reject updating variables which are not
                    % simulated (e.g. "d", "z")
                end

                % use the k-th step references
                for i = 1:length(references)
                    varargin{references(i).position} = references(i).value(:, k);
                end

                % -------------------  initialize MPC  --------------------
                % Sim Step 1
                if k == 1
                    if obj.controller.CRparams.options.CRfixed
                         obj.controller.InActiveCons = obj.controller.CRparams.optprob.nineq;
                    else
                        obj.controller.InitializeCRMPC(obj.controller.CRparams.options.initializeType);
                    end
                end
                % ----> TODO: Reference/uprev changes => reinitialize MPC
                % ---------------------------------------------------------

                [u, feasible, openloop] = obj.controller.evaluate(x0, varargin{:});
                if ~feasible
                    warning('No control action found at step %d, aborting.', k);
                    break
                end

                if isempty(pred_ctrl)
                    pred_ctrl = openloop;
                else
                    pred_ctrl = [pred_ctrl, openloop];
                end


                % note: we could use obj.system.update(u) here, but that
                % introduces significant overhead. update_equation is much
                % faster.
                if include_disturbance
                    [x0, y, ~, d] = obj.system.update_equation(x0, u);
                else
                    [x0, y] = obj.system.update_equation(x0, u);
                end
                X = [X x0];
                U = [U u];
                Y = [Y y];
                if include_disturbance
                    D = [D d];
                end
                J = [J openloop.cost];

                SolverTime(k) = openloop.solvertime;                        % updated
                ncons(k) = openloop.ncons;                                  % updated
                Zdata(:,k) = openloop.Z;                                    % updated
            end
            %--------------------------------------------------------------
            out.X = X;
            out.U = U;
            out.Y = Y;
            if include_disturbance
                out.D = D;
            end
            out.cost = J;

            % Assign further outputs
            out.predictions.controller = pred_ctrl;                         % updated
            out.Z = Zdata;                                                  % updated
            out.ncons = ncons;                                              % updated
            out.solvertime = SolverTime;                                    % updated

            if exist('mpt.extras.DataRecorder', 'file')
                out = mpt.extras.DataRecorder(out);
            end
        end
    end
end