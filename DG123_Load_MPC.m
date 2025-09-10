%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 文件名称: DG123_Load_MPC.m
% 作    用: 为 IEEE123 节点系统的日内动态校正 (MPC) 第一阶段
%           (stage2_1_WithVSM_IEEE123.m) 提供:
%             - 光伏与负荷 96×96 (t,b) 形式的滚动预测矩阵
%             - 节点装机比例 / 负荷比例 / 功率因数角
%             - 日前 OLTC 档位与电容器投切方案并拓展到日内 (含多出的缓冲时段)
%             - 其它辅助基线数据
%
% 升级说明(基于原 DG123_Load.m):
%   1) 保留原 24 点日曲线构造与标幺缩放逻辑
%   2) 生成 96 点 (15min) 基线 (Solar / Load)
%   3) 构造滚动预测矩阵 Solar_data_15, Load_data_15 (96×96)
%   4) 引入可调随机扰动 (可关闭)
%   5) 加入用户提供的“日前值 → 日内展开”代码 (Ktij / NtCB)，并保持可自行修改
%      注意：用户方案把 24 点扩展到 96 点后再补 16 点 (得到 112) 作为 MPC 窗口的
%            安全缓冲 (例如最后几个时刻仍可取到未来预测)。本文件保留该逻辑。
%
% 关键输出变量 (供主调度脚本使用):
%   Solar_radio (123×1)
%   Load_radio  (123×1)
%   q_Load_radio(123×1)
%   p_load, q_load
%   Solar_data (1×24)        — 24点基线
%   Load_data  (1×24)
%   q_Load_data(1×24)
%   Solar_base_96 (1×96)     — 96点基线
%   Load_base_96  (1×96)
%   q_Load_base_96(1×96)
%   Solar_data_15 (96×96)    — 滚动预测矩阵: 行=t当前, 列=未来绝对时刻
%   Load_data_15  (96×96)
%   Wind_data_15  (96×96 这里为零)
%   Ktij_1 (1×112)           — 日内档位序列(含缓冲)
%   NtCB_1, NtCB_2, NtCB_3 (1×112)
%
% 若不需要缓冲，可自行裁剪前 96 列:
%   Ktij_1_use = Ktij_1(1:96);
%
% 使用建议:
%   - 如需关闭随机性: 将 forecast_error_solar / load 和 error_growth_factor_* 设为 0
%   - 可在本脚本最末尾再添加自定义情景扩展
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ================== 基础数据 (需包含 Bus / Branch) ====================
run system123_data.m;  % Bus: 假设列2(kW)列3(kVar); Branch: 拓扑

%% ================== 光伏原始 24 点日曲线 ================================
Solar_origin_data=[113.778,113.778,113.778,113.778,113.778,113.778,...
                   113.778,123.667,157.333,233.333,320.556,366.778,...
                   375.778,366.667,320.667,260.111,160.111,133.222,...
                   113.778,113.778,113.778,113.778,113.778,113.778];

% 光伏装机比例 (与原脚本一致; 基准总装机 7000 作为参考)
Solar_radio=zeros(123,1);
Solar_radio(28)=800/7000;
Solar_radio(33)=900/7000;
Solar_radio(42)=800/7000;
Solar_radio(86)=900/7000;
Solar_radio(92)=600/7000;
Solar_radio(97)=700/7000;
Solar_radio(108)=800/7000;
Solar_radio(111)=700/7000;
Solar_radio(116)=800/7000;

%% ================== 负荷原始 24 点日曲线 ================================
Load_origin_data=[139.778 135.333 129.222 133.222 149.222 197.778,...
                  224.556 204.889 193.556 173.889 174.667 173.556,...
                  169.222 168.667 169.222 183.556 196.111 219.333,...
                  224.778 239.778 226.222 190.667 165 131];

%% ================== 节点基值 (MW) 与比例 =================================
p_load = Bus(:,2)/1000;   % kW -> MW
q_load = Bus(:,3)/1000;   % kVar -> MVar

Load_radio    = zeros(123,1);
q_Load_radio  = zeros(123,1);
for a=2:123
    Load_radio(a)   = p_load(a)/sum(p_load);
    q_Load_radio(a) = q_load(a)/sum(q_load);
end

%% ================== 缩放公式 (与原脚本保持) =============================
alpha = 2.92;
% 光伏基线(24点)
Solar_data = (Solar_origin_data-113.778)/(765.778-113.778) * 4 * alpha * 1.45 * 7/7;
% 负荷基线(24点)
Load_data  = (Load_origin_data-57.2222)/(241.2222-57.2222)*5.6*1;
% 无功 24点
q_Load_data = Load_data*(2710/4885);

%% ================== 24 → 25 (闭合) 用于插值 =============================
time_25     = 1:25;
Solar_25    = [Solar_data Solar_data(1)];
Load_25     = [Load_data Load_data(1)];
q_Load_25   = [q_Load_data q_Load_data(1)];

%% ================== 生成 96 点 (15min) 基线 =============================
% 96 个 15min 时段对应 24 小时
time_96 = linspace(1,25,96+1); % 生成 97 个点(含25)
time_96 = time_96(1:96);       % 取前 96
Solar_base_96  = interp1(time_25, Solar_25,  time_96, 'pchip');
Load_base_96   = interp1(time_25, Load_25,   time_96, 'pchip');
q_Load_base_96 = interp1(time_25, q_Load_25, time_96, 'pchip');

% 节点级基线 (NBUS×96)
p_Solar_base_96      = Solar_radio * Solar_base_96;
p_Load_base_96       = Load_radio  * Load_base_96;
q_Load_base_full_96  = q_Load_radio * q_Load_base_96;

%% ================== 功率因数角 (静态) ===================================
theta = zeros(123,1);
for j=1:123
    if p_load(j)==0
        theta(j)=0;
    else
        theta(j)=atan(q_load(j)/p_load(j));
    end
end
theta_node = repmat(theta',1,24); % 保留原 24 点角度矩阵 (兼容性)

%% ================== 构造滚动预测矩阵 (96×96) ===========================
% Solar_data_15(t,b): 在真实时刻 t 对绝对时段 b 的光伏预测
% Load_data_15(t,b) : 同理
% 引入随机误差 (可设置为0 关闭)
rng(1);  % 固定种子确保复现

forecast_error_solar        = 0.02;   % 基准相对误差
forecast_error_load         = 0.015;
error_growth_factor_solar   = 0.003;  % 随前瞻步长递增
error_growth_factor_load    = 0.002;

Solar_data_15 = zeros(96,96);
Load_data_15  = zeros(96,96);
Wind_data_15  = zeros(96,96);  % 本案例默认无风电

for t_now = 1:96
    for b = 1:96
        k_ahead = max(b - t_now,0);
        sigma_s = forecast_error_solar + k_ahead*error_growth_factor_solar;
        sigma_l = forecast_error_load  + k_ahead*error_growth_factor_load;
        
        base_s = Solar_base_96(b);
        base_l = Load_base_96(b);
        
        err_s = sigma_s * base_s * randn;
        err_l = sigma_l * base_l * randn;
        
        Solar_data_15(t_now,b) = max(base_s + err_s, 0);
        Load_data_15(t_now,b)  = max(base_l + err_l, 0);
        Wind_data_15(t_now,b)  = 0;
    end
end

% （可选）若希望预测无偏，可在此对每行进行缩放校正，这里略。

%% ================== 用户指定的“日前值 → 日内值” (OLTC / 电容) ==========
% 保留用户原始脚本方式, 并说明：
% 日前 24 点档位/投切：Ktij (1×24), NtCB(3×24)
% 通过 repmat 每小时复制 4 次 → 96 点，再额外补 16 点 (→ 112) 作为滚动窗口缓冲
%
% 若主调度 N=8, 外层 t=96 时需要访问到 t+N-1=103, 故长度 ≥103 即可, 112 足够。
% 如果修改 N 或滚动长度，请相应增减缓冲部分。
%
% ---------------- 日前值 ----------------
Ktij=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];      % 1×24
NtCB=[ 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3;       % 3×24
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

% ---------------- 日内展开 (每小时4段) ----------------
Ktij = reshape(repmat(Ktij, 4, 1), 1, []);        % 1×96
Ktij = [Ktij, Ktij(:,1:16)];                      % 1×112 (补16段缓冲)

NtCB = reshape(repmat(NtCB, 4, 1), 3, []);        % 3×96
NtCB = [NtCB, NtCB(:,1:16)];                      % 3×112

% 命名与主调度脚本期望一致
Ktij_1 = Ktij;                % 1×112
NtCB_1 = NtCB(1,:);           % 1×112
NtCB_2 = NtCB(2,:);           % 1×112
NtCB_3 = NtCB(3,:);           % 1×112

% 如果只想用前 96 段:
% Ktij_1_use  = Ktij_1(1:96);
% NtCB_1_use  = NtCB_1(1:96);
% NtCB_2_use  = NtCB_2(1:96);
% NtCB_3_use  = NtCB_3(1:96);

%% ================== 旧 24 点形式下的节点功率(保持可用) ================
p_Solar_24 = zeros(123,24);
p_Load_pre = zeros(123,24);
q_Load_pre = zeros(123,24);
for a=1:24
    p_Solar_24(:,a) = Solar_radio * Solar_data(a);
    p_Load_pre(:,a) = Load_radio  * Load_data(a);
    q_Load_pre(:,a) = q_Load_radio * q_Load_data(a);
end

%% ================== 可选绘图 (按需启用) ================================
% figure; plot(Solar_data,'r-*','LineWidth',2); hold on
% plot(Load_data,'g-*','LineWidth',2);
% legend('光伏(24h)','负荷(24h)'); xlabel('时间/h'); ylabel('有功功率/MW');
%
% figure; plot(Solar_base_96,'r-','LineWidth',1.5); hold on
% plot(Load_base_96,'b-','LineWidth',1.5);
% legend('光伏基线96','负荷基线96'); xlabel('15min 时段'); ylabel('功率/MW');
%
% figure; imagesc(Solar_data_15); colorbar; title('Solar\_data\_15 (t,b)');
% xlabel('b (未来绝对时段)'); ylabel('t (当前时段)');
%
% figure; imagesc(Load_data_15); colorbar; title('Load\_data\_15 (t,b)');

%% ================== 输出变量说明 =======================================
% 供 stage2_1_WithVSM_IEEE123.m 使用的关键变量已经在工作区:
%   Solar_radio, Load_radio, q_Load_radio, p_load, q_load
%   Solar_data_15, Load_data_15, Wind_data_15
%   Ktij_1, NtCB_1, NtCB_2, NtCB_3
%   (基线: Solar_data, Load_data, q_Load_data, 以及 96 点基线)
%
% 如需保存独立数据文件，可自行添加:
% save DG123_Load_MPC_prepared.mat Solar_radio Load_radio q_Load_radio ...
%      p_load q_load Solar_data Load_data q_Load_data ...
%      Solar_base_96 Load_base_96 q_Load_base_96 ...
%      Solar_data_15 Load_data_15 Wind_data_15 ...
%      Ktij_1 NtCB_1 NtCB_2 NtCB_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%