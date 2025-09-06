%% 光伏数据
Solar_origin_data=[113.778,113.778,113.778,113.778,113.778,113.778,...
                    113.778,123.667,157.333,233.333,320.556,366.778,...
                    375.778,366.667,320.667,260.111,160.111,133.222,...
                    113.778,113.778,113.778,113.778,113.778,113.778];
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
%% 负荷数据 57.2222 241.222
Load_origin_data=[139.778 135.333 129.222 133.222 149.222 197.778,...
                  224.556 204.889 193.556,173.889 174.667 173.556,...
                  169.222 168.667 169.222 183.556 196.111 219.333,...
                  224.778 239.778 226.222 190.667 165 131];
Load_radio=zeros(123,1);
q_Load_radio=zeros(123,1);
p_load=Bus(:,2)/1000;
q_load=Bus(:,3)/1000;
for a=2:123
    Load_radio(a)=p_load(a)/sum(p_load);
    q_Load_radio(a)=q_load(a)/sum(q_load);
end
alpha=2.92;
Solar_data=(Solar_origin_data-113.778)/(765.778-113.778)*4*alpha*1.45*7/7;
Load_data=(Load_origin_data-57.2222)/(241.222-57.2222)*5.6*1;
q_Load_data=Load_data*(2710/4885);
%% 插值
MPC_T1=0:1:24-1;
Solar_MPC_data=interp1(1:25,[Solar_data Solar_data(1)],MPC_T1,'pchip');
Load_MPC_data=interp1(1:25,[Load_data Load_data(1)],MPC_T1,'pchip');
q_Load_MPC_data=Load_MPC_data*(2710/4885);
for j=1:123
    if(p_load(j)==0)
        theta(j)=0;
    else
        theta(j)=atan(q_load(j)./p_load(j));
    end
end
theta_node=repmat(theta',1,24);
for a=1:24
    p_Solar(:,a)=Solar_radio*Solar_data(a);
    p_Load_pre(:,a)=Load_radio*Load_data(a);
    q_Load_pre(:,a)=q_Load_radio*q_Load_data(a);
end
plot(Solar_data,'r-*','LineWidth',2);
hold on
plot(Load_data,'g-*','LineWidth',2);
legend('光伏','负荷');
xlabel('时间/h');
ylabel('有功功率/MW');

