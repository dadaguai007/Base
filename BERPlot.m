classdef BERPlot < handle
    % copyright Tianwai@OCG,BIT
    % revision history
    %   2021/11/04: create
    %   2021/11/05: add multiplot functions; add check if rop is increasing
    %   2021/12/29: add support for display minor ticks in the figure;
    %   2022/01/13: add support for adding legend
    % TODO：
    % 
    
    properties
        % configuration of figures
        Config; % 存储图形的配置参数
        flagWhiteMarkerFace = 0; %控制标记的填充颜色
        % interpolation setting
        flagInterpolation = 1;% 标志是否进行插值
        orderInterpolation = 1; 
        % auto zoom configuration
        flagAutoZoom = 1;% 标志是否自动缩放图像。
        % minor tick
        flagMinorTick = 1;% 标志是否显示次要刻度
        % legend
        flagAddLegend = 1;% 标志是否添加图例
    end

    properties (Constant = true)
        BER_axis = [ 0.5 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10...
            1e-11 1e-12];
        color = 'rbkmqcy';
        marker = 'so^d+v*p';
    end


    methods
        function obj = BERPlot()
            % initialize the figure configuration
            Config.AxisLineWidth = 1.25;
            Config.AxisFontSize = 8;
            Config.LineWidth =  1.5;
            Config.MarkerSize = 8;
            Config.FontSize = 12;
            Config.FigureSize = [420 350];
            %Config.AxisLineWidth：坐标轴线的宽度。
            %Config.AxisFontSize：坐标轴上标签的字体大小。
            %Config.LineWidth：曲线的线宽。
            % Config.MarkerSize：标记的大小。
            % Config.FontSize：图像中的字体大小。
            % Config.FigureSize：绘图窗口的大小，以像素为单位。
            obj.Config = Config;
        end
        
        function h = plot(obj,rop,ber,colorId,markerID)
            % Single Plot
            % rop: vector for received optical power
            % ber: vector for bit error rate
            if nargin < 4
                colorId = 1;
                markerID = 1;
            end
            % make rop and ber in row factor
            rop = reshape(rop,1,[]);
            ber = reshape(ber,1,[]);
            % check if the rop vector is in monotonically increasing order
            %检查 rop 向量是否是递增的，如果不是，会将 rop 和 ber 向量进行反转，并显示警告信息。
            if any(diff(rop)<0)
                rop = fliplr(rop);
                ber = fliplr(ber);
                warning('Input rop vector should be in increasing order!');
            end
            % calculate Q-factor from BER
            qfactor_db = obj.ber2q(ber);
            % calculate the range of Q axis
            q_axis_db = obj.ber2q(obj.BER_axis);
            % find the min. and max. for x axis
            x_min = floor(min(rop));
            x_max = ceil(max(rop));
            % find the min. and max. for y axis
            y_min = max(q_axis_db(q_axis_db<min(qfactor_db)));
            y_max = min(q_axis_db(q_axis_db>max(qfactor_db)));

            % now plot in Q factor
            if obj.flagWhiteMarkerFace
                mfc = 'w';
            else
                mfc = obj.color(colorId); 
            end
            
            figure(99);
            h(1) = plot(rop,qfactor_db,...
                    [obj.color(colorId),obj.marker(markerID)],...
                    'Linewidth',obj.Config.LineWidth,...
                    'MarkerEdgeColor',obj.color(colorId),...
                    'MarkerFaceColor',mfc,...
                    'MarkerSize',obj.Config.MarkerSize);
            
            % interpolation if needed
            if obj.flagInterpolation
                figure(99); hold on;
                
                fit_axis = x_min-5:x_max+5;
                p = polyfit(rop,qfactor_db,obj.orderInterpolation);
                q = polyval(p,fit_axis);
                h(2) = plot(fit_axis,q,...
                    ['-',obj.color(colorId)],...
                    'Linewidth',obj.Config.LineWidth);
                % plot a new invisible line only for generating legend
                h(3) = plot(fit_axis+50,q,...
                    ['-',obj.color(colorId),obj.marker(markerID)],...
                    'Color',obj.color(colorId),...
                    'Linewidth',obj.Config.LineWidth,...
                    'MarkerEdgeColor',obj.color(colorId),...
                    'MarkerFaceColor',mfc,...
                    'MarkerSize',obj.Config.MarkerSize);
            end
            % set the axis
            set(gca,'xlim',[x_min-1, x_max+1],...
                'xtick', x_min-1:x_max+1);
            set(gca,'YDir', 'reverse');
            set(gca,'ylim',[y_min y_max],...
                'linewidth' , 1.5,...
                'ytick',q_axis_db,...
                'yticklabel',...
                {'10^0','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}',...
                '10^{-6}','10^{-7}','10^{-8}','10^{-9}','10^{-10}',...
                '10^{-11}','10^{-12}'},...
                'fontsize',obj.Config.FontSize)
            if obj.flagMinorTick
                obj.addYMinorTick;
                set(gca,'YMinorGrid','on','MinorGridAlpha',0.1);
            end
            set(gcf,'pos',[200,200,obj.Config.FigureSize]);
            box on;
            grid on;
            xlabel('Received optical power (dBm)');
            ylabel('Bit error ratio');
            hold off;
        end
        
        function h = multiplot(obj,ropMat,berMat,legendStrArry)
            % Plot multiple BER curves in one figure
            % ropMat: N x M matrix, where N is number of ROP vectors
            % berMat: N x M matrix, where N is number of ROP vectors
            % legendStrArry: an array of string for legend
            if nargin < 4
                legendStrArry = [];
            elseif length(legendStrArry) ~= size(berMat,1)
                error('Legend string should have the same size of berMat!');
            end
            
            if size(ropMat,2) == size(berMat,2) && size(ropMat,1) == 1
                ropMat = repmat(ropMat,size(berMat,1),1);
            end
            
            if ~isequal(size(ropMat), size(berMat))
                error('ropMat and berMat should have the same size!');
            end
            
            [nPlots,~] = size(ropMat);
            
            for idx = 1:nPlots
                % calculate the colorID and markerID
                colorIDvec = 1:ceil(nPlots/length(obj.color))*length(obj.color);
                colorIDvec = colorIDvec(1:nPlots);
                
                markerIDvec = 1:ceil(nPlots/length(obj.marker))*length(obj.marker);
                markerIDvec = markerIDvec(1:nPlots);
                % plot
                figure(99);hold on;
                h{idx} = obj.plot(ropMat(idx,:),berMat(idx,:),colorIDvec(idx),...
                    markerIDvec(idx));
                hold off;
            end
            
            if obj.flagAddLegend
                if isempty(legendStrArry)
                    fprintf('No input legend string, default will be used!\n');
                    % TODO
                end
                for idx = 1:length(h)
                    lv(idx) = h{idx}(3);
                end
               legend(lv,legendStrArry{:} );
            end
        end
        
        function addYMinorTick(obj)
           % generate the minor ticks
           minorTickBER = obj.BER_axis(3:end).'*(9:-1:2);
           minorTickBER = [(4:-2).'*obj.BER_axis(2);...
               reshape(minorTickBER.',[],1)];
           minorTickQ = obj.ber2q(minorTickBER);
           figure(99);
           hA = get(gca);
           hA.YAxis.MinorTickValues = minorTickQ;
           hA.YAxis.MinorTick = 'on';
        end
        
        function q = ber2q(obj,BER) %将BER转换为Q因子（估计信号质量的度量）
            % method 1: Taylor expansion
            q = 20*log10(0.71839 - 1.00208*log10(BER) - 0.08037*log10(BER).^2 - ...
                0.0054*log10(BER).^3 - 2.07568e-4*log10(BER).^4 - 3.36139e-6*log10(BER).^5);
            % method 2: built-in function
            %   this method cannot handle the case: BER = 0.5;
            % q = pow2db(erfcinv(BER.*2)*sqrt(2))*2
        end
    end


end