classdef TubeMPCControllerPlotes < handle
    % Class for creating plots

    methods

        function PlotStateTrajectory(obj,TMPC,x0,Nsim)
            % -------------------------------------------------------------
            %               PlotStateTrajectory(TMPC,x0,Nsim)
            % -------------------------------------------------------------
            % Function plots open-loop / closed-loop state trajectory.
            % Closed-loop is chosen if Nsim is provided.
            % -------------------------------------------------------------

            % First check if we can draw the graph
            if TMPC.model.nx > 3
                error('MPTplus: Can not plot higher dimensional problem than 3!')
            end

            % plot open-loop (closed-loop if Nsim was provided)
            if nargin < 4
                obj.PlotOpenLoop(TMPC,x0)
            else
                obj.PlotClosedLoop(TMPC,x0,Nsim)
            end
        end
    end



    methods (Access = private)

        function PlotOpenLoop(obj,TMPC,x0)

            % Evaluate the controller
            [~,feas,data]  = TMPC.evaluate(x0);
            if ~feas
                error('MPTplus: Infeasible initial condition was provided!')
            end

            figure
            hold on
            box on
            grid on

            % -------- Plot nominal trajectory (axis initialization) ----------
            obj.plotState(TMPC,data.X,'b-*',2)

            % -------- Plot terminal set ----------
            if TMPC.model.nx == 1
                Xf = TMPC.TMPCparams.params_TVTMPC.Xf;
                xmin = Xf.Internal.lb;
                xmax = Xf.Internal.ub;
                yline(xmin,'Color',[0.7 0.7 0.7],'LineWidth',2)
                yline(xmax,'Color',[0.7 0.7 0.7],'LineWidth',2)
            else
                TMPC.TMPCparams.params_TVTMPC.Xf.plot('Color',[0.7, 0.7, 0.7])
            end

            % -------- Plot time-varying tubes ----------
            for k = 2:TMPC.TMPCparams.N+1
                if k < TMPC.TMPCparams.Ns
                    Tube = TMPC.ConstructExplicitTube(k-1);
                end
                if TMPC.model.nx == 1
                    % 1D: Tube is an interval [lb, ub], plot as vertical line or band
                    xmin = min(Tube.V) + data.X(:,k);
                    xmax = max(Tube.V) + data.X(:,k);
                    % Option 1: vertical lines
                    plot([k k],[xmin xmax],'Color','red','LineWidth',1.5)
                    % Option 2 (better): filled patch
                    patch([k-0.4 k+0.4 k+0.4 k-0.4],[xmin xmin xmax xmax], ...
                        [1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.3)

                elseif TMPC.model.nx == 2
                    % 2D: Tube + trajectory is a polygon, plot normally
                    plot(Tube + data.X(:,k),'Color','red')
                else
                    % 3D: optional, plot as polygon in 3D space
                    plot3(Tube(1,:) + data.X(1,k), Tube(2,:) + data.X(2,k), Tube(3,:) + data.X(3,k), 'Color','red')
                end
            end

            % -------- Plot trajectory again ----------
            obj.plotState(TMPC,data.X,'b-*',1)

            obj.labelStateAxes(TMPC)
            title('Nominal state trajectory')
            legend({'$\hat{x}_0$','$X_f$','$\hat{T}_k$'}, ...
       'Interpreter','latex','Location','best')
        end



        function PlotClosedLoop(obj,TMPC,x0,Nsim)

            data = TMPC.simulate(x0,Nsim,zeros(TMPC.model.nx,Nsim));

            figure
            hold on
            box on
            grid on

            if TMPC.isExplicit
                TMPC.partition.Domain.plot('Color','yellow')
            end

            % -------- Plot nominal and real trajectory ----------
            obj.plotState(TMPC,data.Xnominal,'b--*',2)
            obj.plotState(TMPC,data.X,'k--*',1)

             % -------- Plot terminal set ----------
            if TMPC.model.nx == 1
                Xf = TMPC.TMPCparams.params_TVTMPC.Xf;
                xmin = Xf.Internal.lb;
                xmax = Xf.Internal.ub;
                yline(xmin,'Color',[0.7 0.7 0.7],'LineWidth',2)
                yline(xmax,'Color',[0.7 0.7 0.7],'LineWidth',2)
            else
                TMPC.TMPCparams.params_TVTMPC.Xf.plot('Color',[0.7, 0.7, 0.7])
            end

            % -------- Plot shifted terminal sets ----------
            if TMPC.model.nx == 1
                % 1D: plot shifted terminal set as horizontal band at each time step
                Xf = TMPC.TMPCparams.params_TVTMPC.Xf;
                xmin = min(Xf.V);
                xmax = max(Xf.V);

                for k = 1:Nsim
                    xnom = data.Xnominal(1,k);
                    % Option 1: vertical lines
                    plot([k k],[xmin+xnom xmax+xnom],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
                    % Option 2 (better): filled patch
                    patch([k-0.4 k+0.4 k+0.4 k-0.4],[xmin+xnom xmin+xnom xmax+xnom xmax+xnom], ...
                        [0.85 0.85 0.85],'EdgeColor','none','FaceAlpha',0.3)
                end
            else
                % 2D / 3D: plot normally
                for k = 1:Nsim
                    plot(TMPC.TMPCparams.params_TVTMPC.Xf + data.Xnominal(:,k), 'Color', [0.7, 0.7, 0.7])
                end
            end

            % -------- Plot again (on top) ----------
            obj.plotState(TMPC,data.Xnominal,'b--*',1)
            obj.plotState(TMPC,data.X,'k--*',1)

            obj.labelStateAxes(TMPC)
            title('Closed-loop state trajectory')

            if TMPC.isExplicit
                legend({'$C_{\infty}$','$\hat{x}_0$','$x(t)$','$X_\text{f}$'}, ...
       'Interpreter','latex','Location','best')
            else
                legend({'$\hat{x}_0$','$x(t)$','$X_f$'}, ...
       'Interpreter','latex','Location','best')
            end
        end



        function plotState(obj,TMPC,X,style,width)
            switch TMPC.model.nx
                case 1
                    plot(X(1,:),style,'LineWidth',width)
                case 2 
                   plot(X(1,:),X(2,:),style,'LineWidth',width)
                case 3
                    plot3(X(1,:),X(2,:),X(3,:),style,'LineWidth',width)
            end
        end



        function labelStateAxes(obj,TMPC)
            switch TMPC.model.nx
                case 1
                    xlabel('time step','Interpreter','latex')
                    ylabel('$x_1$','Interpreter','latex')
                case 2
                    xlabel('$x_1$','Interpreter','latex')
                    ylabel('$x_2$','Interpreter','latex')
                case 3
                    xlabel('$x_1$','Interpreter','latex')
                    ylabel('$x_2$','Interpreter','latex')
                    zlabel('$x_3$','Interpreter','latex')
                    view(3)
            end
        end

    end
end