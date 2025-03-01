% predistortion have the digital and analog way to implement
% the function is only to the mzm
% Vbias = 0.5; % bias voltage normalized by Vpi
% Vswing = 1;  % normalized voltage swing. 1 means that modulator is driven at full scale

% 进行预失真前，还要先确定原先的电平值和判决值即：throad 和 Vk
function self = mzm_predistortion(M, Vswing, Vbias,type,xd)
L = reset_levels(M);
self=unbias(L);
throad=self.b;
switch(lower(type))
    case 'digital'
        %% Predistort levels to compensate for MZM nonlinear response in IM-DD
        % get the bias V from the power
        predist = @(p) 2/pi*asin(sqrt(p));
        % nonlinear response
        dist = @(v) sin(pi/2*v)^2;

        Vmin = Vbias - Vswing/2;
        Vmax = Vbias + Vswing/2;
        %这里，dist 函数用于计算输入电压对应的失真功率。
        % Pmax 和 Pmin 是在输入电压范围内得到的失真功率的上限和下限。
        Pmax = dist(Vmax);
        Pmin = dist(Vmin);
        %输入电压范围内，对应于失真功率范围的步进
        DP = (Pmax-Pmin)/(M-1);
        %对应了预失真的功率值
        Pk = Pmin:DP:Pmax;

        % 首先通过失真函数 dist 将输入电压转换为失真功率。
        % 然后，应用预失真函数 predist 对这些失真功率进行反操作，以获得相应的预失真电压。

        % Predistortion
        Vk = predist(Pk);

        % Set
        self = set_levels(Vk, throad);


        % 注:判决阈值没有预失真，因为预失真只用于在发射机产生电平。
        % 最好将没预失真前信号电平图绘制出来
        figure(233), clf, hold on, box on
        t = linspace(0, 1);
        plot(t, sin(pi/2*t).^2, 'k');
        plot((self.a*[1 1]).', [zeros(1, M); sin(pi/2*self.a.').^2], 'k');
        % plot([zeros(1, self.M); self.a.'], ((Pk.')*[1 1]).', 'k')
        xlabel('Driving signal')
        ylabel('Resulting power levels')
        axis([0 1 0 1])

        % 得到 self 之后需要进行 mod signal，原先使用的pammod可能会失效


    case 'analog'
        self.out = 2/pi*asin(sqrt(abs(xd))).*sign(xd); % apply predistortion
end


%% Levels and decision thresholds
    function self = set_levels(levels, thresholds)
        %设置电平（levels）和决策阈值（thresholds）
        %% Set levels to desired values
        %确保输入参数 levels 和 thresholds 是列向量（M x 1），而不是行向量（1 x M）。这是通过对矢量进行转置来实现的
        if size(levels, 1) < size(levels, 2) % ensure levels are M x 1 vector
            levels = levels.';
        end
        if size(thresholds, 1) < size(thresholds, 2)  % ensure thresholds are M x 1 vector
            thresholds = thresholds.';
        end

        self.a = levels;
        self.b = thresholds;
    end

    function self = reset_levels(M)
        %归一化后的取值
        % v_value
        self.a = ((0:2:2*(M-1))/(2*(M-1))).';
        %throad_value
        self.b = ((1:2:(2*(M-1)-1))/(2*(M-1))).';
    end
  % pam信号的初始电平应该是-3——3之间，所以有必要对信号去除直流和归一化，将bit信号转换pam信号
    function self = unbias(self)
        %% Remove DC bias from levels and normalize to have excusion from -1 to 1
        self.b = self.b - mean(self.a);
        self.a = self.a - mean(self.a);
        self.b = self.b/self.a(end);
        self.a = self.a/self.a(end);
    end

    
end


