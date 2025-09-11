clear all
clc

run ieee_33_node_system.m;
run DG_Load_MPC;
run oltc_cbs_and_withoutvsm.m;
tic
t1 = clock;
%% 常量定义
%------------------ 基准值设置--------------------------------------
U_b=12.66;%kV
S_b=1;%MW
Z_b=U_b^2/S_b;
%----------------电阻电抗--------------------------------
r_ij = Branch( : , 4 )/Z_b;
x_ij = Branch( : , 5 )/Z_b;
%------------------目标函数值--------------------------------------
N =8;%时间
result = zeros(1,N);
%------------------电压上下限--------------------------------------
a1=0.95;
a2=1.05;
%----------------变压器抽头及电容器--------------------------------
delta_T = 1;%hour
K_max = 5;%抽头最大挡位
kk = (0:K_max*2);
k_ij0 = 1;%初始变比

delta_kij = 0.1/(K_max*2);%每个挡位电压变化率
%------------------电容器单位无功功率--------------------------------------
qCB = 0.06;%每个电容器的无功功率（Mvar）

%----------------光伏安装位置及其容量（Mwh）--------------------------------
q_n=6;

pv1=10;
pv2=14;
pv3=17;
pv4=24;
pv5=28;%光伏安装位置
pv6=32;

S_solar=zeros(1,q_n);
S_solar(1)=0.8*1.1;%光伏容量
S_solar(2)=0.3*1.1;
S_solar(3)=0.3*1.1;
S_solar(4)=2.3*1.1;
S_solar(5)=3.5*1.1;
S_solar(6)=1*1.1;
%----------------光伏无功出力限制值--------------------------------
Q_pv=zeros(1,q_n);
for q_pv=1:q_n
    Q_pv(q_pv)=S_solar(q_pv)*0.436;
end
%----------------SOP无功出力限制值--------------------------------
Q_sop1 = zeros(4,112);

%   for  j=1:112
%    Q_sop1(1 , j)= sqrt(4-P1_data20(1 , j)^2);
%    Q_sop1(2 , j)= sqrt(4-P1_data20(2 , j)^2);
%    Q_sop1(3 , j)= sqrt(4-P2_data20(1 , j)^2);
%    Q_sop1(4 , j)= sqrt(4-P2_data20(2 , j)^2);
%   end
%  
% save sop无功出力最大值_112 Q_sop1; 
%---------------------ESSOP各子系统容量和损耗系数-------------------------------
%------容量------------
S_acdc = 2;
S_dcdc = 0.6;
S_ess = 1;
%------损耗系数------------
C_a0=0.005;
C_a1=0.005;
C_a2=0.02;

C_d0=0.0035;
C_d1=0.0025;
C_d2=0.014;

% C_aux=0.002*S_ess;
C_aux=0;
C_ES=0.016;
%-----------------------SOP有功无功、PV无功参考值、开关动作-----------------------------------    
% 96组数据，依次代入进行求解
load('日内步骤1传递给2的数据1_1_112_WithoutVSM_随机性.mat');
load('风光负荷_96_112_15min_随机性.mat');
load('sop无功出力最大值_112.mat');
%% 多时间尺度循环
for t=1:96
 b1 = zeros(2,N); 
 b2 = zeros(2,N);

 P1=zeros(2,N);
 P2=zeros(2,N);

 Ktij = zeros(1,N);
 NtCB =zeros(3,N);
 
 Qsop_ref=zeros(4,N);
 Qpv_ref=zeros(6,N);

 t_1=t;
 t_2=t+N-1;  

    Ktij( 1 ,:)= Ktij_1( 1,t_1:t_2 ); 
    
    NtCB( 1 ,:)= NtCB_1( 1,t_1:t_2 ); 
    NtCB( 2 ,:)= NtCB_2( 1,t_1:t_2 );
    NtCB( 3 ,:)= NtCB_3( 1,t_1:t_2 );  


    P1( 1 , :)=P1_data20( 1,t_1:t_2);  
    P1( 2 , :)=P1_data20( 2,t_1:t_2);
    P2( 1 , :)=P2_data20( 1,t_1:t_2);
    P2( 2 , :)=P2_data20( 2,t_1:t_2);
    
    b1( 1 , :)=b1_data20(1,t_1:t_2); 
    b1( 2 , :)=b1_data20(2,t_1:t_2);  
    b2( 1 , :)=b2_data20(1,t_1:t_2); 
    b2( 2 , :)=b2_data20(2,t_1:t_2);  
    
   Qsop_ref( 1 , :)=Q1_data20( 1,t_1:t_2); 
   Qsop_ref( 2 , :)=Q1_data20( 2,t_1:t_2);  
   Qsop_ref( 3 , :)=Q2_data20( 1,t_1:t_2);  
   Qsop_ref( 4 , :)=Q2_data20( 2,t_1:t_2);  
   
   Qpv_ref( 1 , :)= q_Solar_data20( 1,t_1:t_2); 
   Qpv_ref( 2 , :)= q_Solar_data20( 2,t_1:t_2); 
   Qpv_ref( 3 , :)= q_Solar_data20( 3,t_1:t_2); 
   Qpv_ref( 4 , :)= q_Solar_data20( 4,t_1:t_2); 
   Qpv_ref( 5 , :)= q_Solar_data20( 5,t_1:t_2); 
   Qpv_ref( 6 , :)= q_Solar_data20( 6,t_1:t_2); 
   
   Q_sop( 1 , :)=Q_sop1( 1,t_1:t_2); 
   Q_sop( 2 , :)=Q_sop1( 2,t_1:t_2); 
   Q_sop( 3 , :)=Q_sop1( 3,t_1:t_2); 
   Q_sop( 4 , :)=Q_sop1( 4,t_1:t_2);

%----------------负荷光伏数据及恒定功率因数--------------------------------
p_Solar = zeros(33,N);
p_Wind  = zeros(33,N);
p_Load  = zeros(33,N);
q_Load  = zeros(33,N);
tan_theta = q_load ./ p_load;%恒定功率因数
for b = t_1:t_2
    p_Solar( : , b-t+1 ) = Solar_radio * Solar_data(t, b );
    p_Load ( : , b-t+1 ) = Load_radio  *  Load_data(t, b );
    q_Load ( : , b-t+1) = ( p_Load ( : , b-t+1 ) ./ p_load ) .* ( q_load );
end

%% 定义初始运行点决策变量
q_Solar =sdpvar(q_n, N, 'full');
%-----------------光伏削减量------------------------------------
% PV_cut=sdpvar(q_n, N, 'full');
PV_cut=zeros(q_n, N);
P_pv=sdpvar(q_n, N, 'full');
%----------------电压、支路有功、支路无功--------------------------------
x_ui_square = sdpvar(33, N, 'full');
x_pij = sdpvar(32, N, 'full');
x_qij = sdpvar(32, N, 'full');
%--------------------1组ESSOP/双馈线-------------------------------------
sv=2;
%AC/DC视在功率辅助变量/非负
k11= sdpvar(sv, N, 'full');
k12= sdpvar(sv, N, 'full');
k21= sdpvar(sv, N, 'full');
k22= sdpvar(sv, N, 'full');
%四个子系统损耗
P1_L=sdpvar(sv, N, 'full');
P2_L=sdpvar(sv, N, 'full');
%AC/DC无功功率
Q1=sdpvar(sv, N, 'full');
Q2=sdpvar(sv, N, 'full');
%-----------------电压偏差------------------------------------
Aux = sdpvar(33, N, 'full');
%----------------损耗-----------------------------
P_loss=sdpvar(32, N, 'full'); 
%% 下垂控制参数
% ----------------光伏-----------------------------
M=1000;
n_i=q_n;
v=sdpvar(n_i, N, 'full');
aa=sdpvar(n_i, N, 6, 'full');
k1=sdpvar(n_i, N, 2, 'full');
dd=binvar(n_i, N, 5, 'full');

bat=sdpvar(n_i, N, 2, 'full');
w1=sdpvar(n_i, N, 4, 'full');
w2=sdpvar(n_i, N, 4, 'full');

L1=binvar(n_i, N, 4, 'full');
L2=binvar(n_i, N, 4, 'full');
%-----------------SOP-----------------------------
M=1000;
S_i=4;
S_v=sdpvar(S_i, N, 'full');
S_aa=sdpvar(S_i, N, 6, 'full');
S_k1=sdpvar(S_i, N, 2, 'full');
S_dd=binvar(S_i, N, 5, 'full');

S_bat=sdpvar(S_i, N, 2, 'full');
S_w1=sdpvar(S_i, N, 4, 'full');
S_w2=sdpvar(S_i, N, 4, 'full');

S_L1=binvar(S_i,  N,4, 'full');
S_L2=binvar(S_i,  N,4, 'full');
%% 定义约束条件
Constraints = [];
% Constraints = [ Constraints , q_Solar== 0 ];
% Constraints = [ Constraints , b1 == 0 ];
% Constraints = [ Constraints , b2 == 0 ];
% Constraints = [ Constraints , b3 == 0 ];
for opt_num = 1: (N)
     
    tic
%% 目标函数

        f(opt_num)=1000*0.25*(0.08*(sum(P_loss( : , opt_num ))+sum( P1_L( : , opt_num )+ P2_L( : , opt_num )))...
                  + 0.64*sum(PV_cut( : , opt_num )))+25*(sum( Aux( : , opt_num )));  
 

%           f(opt_num)=(sum( abs(x_ui_square( : , opt_num )-1)));  
            
        f_start(opt_num)= sum(P_loss( : , opt_num ))*(1/4);    
        f_essop(opt_num)= (sum( P1_L( : , opt_num )+ P2_L( : , opt_num )))*(1/4);
        f_PV_cut(opt_num)= sum(PV_cut( : , opt_num ))*(1/4);
        
%% 平衡节点电压约束
   Constraints = [ Constraints , x_ui_square( 1 , opt_num ) == 1 ];
   
%% 有功无功平衡
for k = 2 : 33
    node_out = find(Branch(:,2) == k);%下游支路
    node_in  = find(Branch(:,3) == k);%上游支路

        Constraints = [ Constraints ,P_loss(k-1 , opt_num ) >= r_ij( k-1 )*(x_pij( k-1 , opt_num )^2 + x_qij( k-1 , opt_num )^2) ]; 
    
    if( k == 2 )
        Constraints = [ Constraints ,x_pij( node_in , opt_num ) + p_Solar( k , opt_num ) + p_Wind( k , opt_num )...
                                     - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , x_qij( 1 , opt_num )...
                                      - q_Load( k, opt_num  )== sum( x_qij( node_out , opt_num ) ) ];
                                  
%% ESSOP1    
    elseif( k == 12 )
        Constraints = [ Constraints , x_pij( node_in , opt_num )  ...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num )...
                                      + P1( 1 , opt_num )- p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , x_qij( node_in ,opt_num )  ...
                                      + Q1( 1 , opt_num )- q_Load( k , opt_num  ) == sum( x_qij( node_out , opt_num ) ) ];
      
    elseif( k == 22 )
        Constraints = [ Constraints , sum( x_pij( node_in ,opt_num ) )...
                                       + P2( 1 , opt_num )+ p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                       + Q2( 1 , opt_num )- q_Load( k , opt_num  )== 0 ];
                     
%% ESSOP2 
    elseif( k == 25 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                     + P1(2 , opt_num )+ p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in ,opt_num )  )...
                                     + Q1( 2 , opt_num )- q_Load( k , opt_num ) == 0 ];                  
                                
    elseif( k == 29 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                     + P2( 2 , opt_num )+ p_Solar( k , opt_num ) + p_Wind( k , opt_num )...
                                    - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                     + Q2( 2 , opt_num )- q_Load( k , opt_num )== sum( x_qij( node_out , opt_num ) ) ];                     
                                
    elseif( k == 18 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                    + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                    - q_Load( k , opt_num ) == 0 ];                        
  
    elseif( k == 33 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                    + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                    - q_Load( k , opt_num  ) == 0 ];

                                                                                              
    elseif( k == 6 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ NtCB(1,opt_num  )*qCB   == sum( x_qij( node_out , opt_num ) ) ];
                                  
    elseif( k == 13 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ NtCB( 2,opt_num   )*qCB   == sum( x_qij( node_out , opt_num ) ) ];
                           

                                                                  
     elseif( k == 31     )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )  + NtCB( 3,opt_num   )*qCB  == sum( x_qij( node_out , opt_num ) ) ];
                                                      
                                         
%% 含风光机组的节点有功无功功率平衡    
    elseif( k == pv1 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                      + p_Solar( k , opt_num )-PV_cut( 1 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ q_Solar(1 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];
                         
                                  
    elseif( k ==pv2 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num )-PV_cut( 2 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ q_Solar(2 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];

    elseif( k ==pv3 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num )-PV_cut( 3 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ q_Solar( 3 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];    
                                                         
   elseif( k ==pv4)
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num )-PV_cut( 4 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ q_Solar(4 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];
                                                                                            
   elseif( k == pv5 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                      + p_Solar( k , opt_num )-PV_cut( 5 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ q_Solar(5 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];
                                  
   elseif( k == pv6 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                      + p_Solar( k , opt_num )-PV_cut( 6 , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ q_Solar( 6 , opt_num ) == sum( x_qij( node_out , opt_num ) ) ];    
                         
    else
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num ) == sum( x_qij( node_out , opt_num ) ) ]; 

    end         
end
  
%% 支路欧姆定律
for r = 2 : 32
    Constraints = [ Constraints , x_ui_square( Branch(r,2) , opt_num ) - x_ui_square( Branch(r,3) , opt_num )  ...
                   -  ( r_ij(r) * x_pij( r , opt_num ) + x_ij(r)* x_qij( r , opt_num ) ) == 0 ];
end
              
%% 运行电压约束
Constraints = [ Constraints , x_ui_square( : , opt_num ) >= a1 ];
Constraints = [ Constraints , x_ui_square( : , opt_num ) <= a2 ];

%% 风光机组无功功率容量和最值约束
for s_n=1:q_n
%------------------------容量约束---------------------------------------
Constraints = [ Constraints ,(q_Solar(s_n , opt_num )^2 + P_pv( s_n , opt_num )^2 )<=2*( S_solar(s_n)/ ( sqrt(2) ) ) * (S_solar(s_n)/ ( sqrt(2) ) ) ];
%------------------------最值约束/功率因数约束---------------------------------------
Constraints = [ Constraints ,-Q_pv(s_n) <= q_Solar(s_n , opt_num )<= Q_pv(s_n)];
end

%------------------------sop功率因数无功约束---------------------------------------
% Constraints = [ Constraints ,-Q_sop <= Q1(: , opt_num )<= Q_sop];
% Constraints = [ Constraints ,-Q_sop <= Q2(: , opt_num )<= Q_sop];

%% 光伏下垂控制
% ------------------------初始运行点---------------------------------------
for q_i=1:6
     if (q_i==1)
Constraints = [ Constraints , q_Solar(q_i , opt_num )== v(q_i , opt_num )*Q_pv(1,q_i)]; 
Constraints = [ Constraints , x_ui_square( 10 , opt_num )== 0.9*aa(q_i , opt_num ,1)+0.92*aa(q_i , opt_num , 2 )+k1(q_i , opt_num ,1)+k1(q_i , opt_num ,2)+1.08*aa(q_i , opt_num , 5 )+1.1*aa(q_i , opt_num , 6 )];
Constraints = [ Constraints , v(q_i , opt_num )== aa(q_i , opt_num , 1)+aa(q_i, opt_num , 2 )+aa(q_i, opt_num , 3 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)+aa(q_i , opt_num , 4 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)...
                               - aa(q_i , opt_num , 5 )- aa(q_i , opt_num , 6 )];  
     elseif( q_i == 2 ) 
Constraints = [ Constraints , q_Solar(q_i , opt_num )== v(q_i , opt_num )*Q_pv(1,q_i)];    
Constraints = [ Constraints , x_ui_square( 14 , opt_num )== 0.9*aa(q_i , opt_num ,1)+0.92*aa(q_i , opt_num , 2 )+k1(q_i , opt_num ,1)+k1(q_i , opt_num ,2)+1.08*aa(q_i , opt_num , 5 )+1.1*aa(q_i , opt_num , 6 )];
Constraints = [ Constraints , v(q_i , opt_num )== aa(q_i , opt_num , 1)+aa(q_i, opt_num , 2 )+aa(q_i, opt_num , 3 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)+aa(q_i , opt_num , 4 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)...
                               - aa(q_i , opt_num , 5 )- aa(q_i , opt_num , 6 )];    
     elseif( q_i == 3 ) 
Constraints = [ Constraints , q_Solar(q_i , opt_num )== v(q_i , opt_num )*Q_pv(1,q_i)];    
Constraints = [ Constraints , x_ui_square( 17 , opt_num )== 0.9*aa(q_i , opt_num ,1)+0.92*aa(q_i , opt_num , 2 )+k1(q_i , opt_num ,1)+k1(q_i , opt_num ,2)+1.08*aa(q_i , opt_num , 5 )+1.1*aa(q_i , opt_num , 6 )];
Constraints = [ Constraints , v(q_i , opt_num )== aa(q_i , opt_num , 1)+aa(q_i, opt_num , 2 )+aa(q_i, opt_num , 3 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)+aa(q_i , opt_num , 4 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)...
                               - aa(q_i , opt_num , 5 )- aa(q_i , opt_num , 6 )];   
     elseif( q_i == 4 ) 
Constraints = [ Constraints , q_Solar(q_i , opt_num )== v(q_i , opt_num )*Q_pv(1,q_i)];    
Constraints = [ Constraints , x_ui_square( 24 , opt_num )== 0.9*aa(q_i , opt_num ,1)+0.92*aa(q_i , opt_num , 2 )+k1(q_i , opt_num ,1)+k1(q_i , opt_num ,2)+1.08*aa(q_i , opt_num , 5 )+1.1*aa(q_i , opt_num , 6 )];
Constraints = [ Constraints , v(q_i , opt_num )== aa(q_i , opt_num , 1)+aa(q_i, opt_num , 2 )+aa(q_i, opt_num , 3 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)+aa(q_i , opt_num , 4 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)...
                               - aa(q_i , opt_num , 5 )- aa(q_i , opt_num , 6 )];  
     elseif( q_i == 5 ) 
Constraints = [ Constraints , q_Solar(q_i , opt_num )== v(q_i , opt_num )*Q_pv(1,q_i)];    
Constraints = [ Constraints , x_ui_square( 28 , opt_num )== 0.9*aa(q_i , opt_num ,1)+0.92*aa(q_i , opt_num , 2 )+k1(q_i , opt_num ,1)+k1(q_i , opt_num ,2)+1.08*aa(q_i , opt_num , 5 )+1.1*aa(q_i , opt_num , 6 )];
Constraints = [ Constraints , v(q_i , opt_num )== aa(q_i , opt_num , 1)+aa(q_i, opt_num , 2 )+aa(q_i, opt_num , 3 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)+aa(q_i , opt_num , 4 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)...
                               - aa(q_i , opt_num , 5 )- aa(q_i , opt_num , 6 )];  
     elseif( q_i == 6 ) 
Constraints = [ Constraints , q_Solar(q_i , opt_num )== v(q_i , opt_num )*Q_pv(1,q_i)];    
Constraints = [ Constraints , x_ui_square( 32 , opt_num )== 0.9*aa(q_i , opt_num ,1)+0.92*aa(q_i , opt_num , 2 )+k1(q_i , opt_num ,1)+k1(q_i , opt_num ,2)+1.08*aa(q_i , opt_num , 5 )+1.1*aa(q_i , opt_num , 6 )];
Constraints = [ Constraints , v(q_i , opt_num )== aa(q_i , opt_num , 1)+aa(q_i, opt_num , 2 )+aa(q_i, opt_num , 3 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)+aa(q_i , opt_num , 4 )*Qpv_ref(q_i , opt_num )/Q_pv(1,q_i)...
                               - aa(q_i , opt_num , 5 )- aa(q_i , opt_num , 6 )];  
     end           
Constraints = [ Constraints , 0<=aa(q_i , opt_num , : )]; 
Constraints = [ Constraints , sum(aa(q_i , opt_num , : ))==1];                          

Constraints = [ Constraints , aa(q_i , opt_num , 1)<=dd(q_i , opt_num , 1)]; 
Constraints = [ Constraints , aa(q_i , opt_num , 6)<=dd(q_i, opt_num , 5)]; 
for d_n=2:5
    Constraints = [ Constraints , aa(q_i, opt_num , d_n)<= dd(q_i , opt_num , d_n)+ dd(q_i , opt_num , d_n-1)]; 
end
Constraints = [ Constraints , sum(dd(q_i , opt_num , : ))==1]; 

Constraints = [ Constraints , k1(q_i , opt_num ,1)==0.92*aa(q_i , opt_num , 3)+0.01*bat(q_i , opt_num ,1)];  
Constraints = [ Constraints , bat(q_i , opt_num ,1)==w1(q_i , opt_num,1 )+2*w1(q_i , opt_num,2 )+4*w1(q_i , opt_num,3 )+8*w1(q_i, opt_num,4 )];
for w_m=1:4
Constraints = [ Constraints , aa(q_i , opt_num , 3)-(1-L1(q_i ,opt_num,w_m ))*M<=w1(q_i , opt_num,w_m )]; 
Constraints = [ Constraints , w1(q_i , opt_num,w_m )<= aa(q_i , opt_num , 3)];

Constraints = [ Constraints , 0<=w1(q_i , opt_num ,w_m )]; 
Constraints = [ Constraints , w1(q_i , opt_num ,w_m )<= L1(q_i ,opt_num,w_m )*M]; 
end

Constraints = [ Constraints , k1(q_i , opt_num ,2)==0.92*aa(q_i , opt_num , 4)+0.01*bat(q_i , opt_num ,2)];  
Constraints = [ Constraints , bat(q_i , opt_num ,2)==w2(q_i , opt_num,1 )+2*w2(q_i , opt_num,2 )+4*w2(q_i , opt_num,3 )+8*w2(q_i, opt_num,4 )];
for w_m=1:4
Constraints = [ Constraints , aa(q_i , opt_num , 4)-(1-L2(q_i ,opt_num,w_m ))*M<=w2(q_i , opt_num,w_m )]; 
Constraints = [ Constraints , w2(q_i , opt_num,w_m )<= aa(q_i , opt_num , 4)]; 
Constraints = [ Constraints , 0<= w2(q_i , opt_num,w_m )]; 
Constraints = [ Constraints , w2(q_i , opt_num,w_m )<= L2(q_i ,opt_num,w_m )*M]; 
end

% Constraints = [ Constraints ,   1<=L1(q_i ,1 )+2*L1(q_i ,2 )+4*L1(q_i ,3 )+8*L1(q_i ,4 )+16*L1(q_i ,5)]; 
% Constraints = [ Constraints ,   L1(q_i ,1 )+2*L1(q_i ,2 )+4*L1(q_i ,3 )+8*L1(q_i,4 )+16*L1(q_i ,5)<=19]; 
% Constraints = [ Constraints ,   1<=L2(q_i,1 )+2*L2(q_i ,2 )+4*L2(q_i ,3 )+8*L2(q_i,4 )+16*L2(q_i ,5)]; 
% Constraints = [ Constraints ,   L2(q_i ,1 )+2*L2(q_i ,2 )+4*L2(q_i ,3 )+8*L2(q_i ,4 )+16*L2(q_i ,5)<=19]; 
% Constraints = [ Constraints ,   (L1(q_i ,1 )+2*L1(q_i ,2 )+4*L1(q_i ,3 )+8*L1(q_i ,4 )+16*L1(q_i ,5))<=(L2(q_i,1 )+2*L2(q_i,2 )+4*L2(q_i,3 )+8*L2(q_i,4 )+16*L2(q_i ,5))]; 
Constraints = [ Constraints ,   1<=L1(q_i ,opt_num,1 )+2*L1(q_i ,opt_num,2 )+4*L1(q_i ,opt_num,3 )+8*L1(q_i ,opt_num,4 )]; 
Constraints = [ Constraints ,   L1(q_i ,opt_num,1 )+2*L1(q_i ,opt_num,2 )+4*L1(q_i ,opt_num,3 )+8*L1(q_i ,opt_num,4 )<=15]; 
Constraints = [ Constraints ,   1<=L2(q_i,opt_num,1 )+2*L2(q_i ,opt_num,2 )+4*L2(q_i ,opt_num,3 )+8*L2(q_i,opt_num,4 )]; 
Constraints = [ Constraints ,   L2(q_i ,opt_num,1 )+2*L2(q_i ,opt_num,2 )+4*L2(q_i ,opt_num,3 )+8*L2(q_i ,opt_num,4 )<=15]; 
Constraints = [ Constraints ,   L1(q_i ,opt_num,1 )+2*L1(q_i ,opt_num,2 )+4*L1(q_i ,opt_num,3 )+8*L1(q_i ,opt_num,4 )<=L2(q_i ,opt_num,1 )+2*L2(q_i ,opt_num,2 )+4*L2(q_i ,opt_num,3 )+8*L2(q_i ,opt_num,4 )]; 

end

%% SOP下垂控制
%------------------------初始运行点---------------------------------------
for s_i=1:4                              
      if( s_i == 1 && b1(1 , opt_num)== 1 )    
         Constraints = [ Constraints , Q1(s_i , opt_num )== S_v(s_i , opt_num )*Q_sop(s_i , opt_num )];    
         Constraints = [ Constraints , x_ui_square( 12 , opt_num )== 0.9*S_aa(s_i , opt_num ,1)+0.92*S_aa(s_i , opt_num , 2 )+S_k1(s_i , opt_num ,1)+S_k1(s_i , opt_num ,2)+1.08*S_aa(s_i , opt_num , 5 )+1.1*S_aa(s_i , opt_num , 6 )];
         Constraints = [ Constraints , S_v(s_i , opt_num )== S_aa(s_i , opt_num , 1)+ S_aa(s_i, opt_num , 2 )+ S_aa(s_i, opt_num , 3 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i , opt_num )+S_aa(s_i , opt_num , 4 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i ,opt_num )...
                                       - S_aa(s_i , opt_num , 5 )- S_aa(s_i , opt_num , 6 )];
                                   
         Constraints = [ Constraints , 0<=S_aa(s_i , opt_num , : )]; 
         Constraints = [ Constraints , sum(S_aa(s_i , opt_num , : ))==1];                          

         Constraints = [ Constraints , S_aa(s_i , opt_num , 1)<=S_dd(s_i , opt_num , 1)]; 
         Constraints = [ Constraints , S_aa(s_i , opt_num , 6)<=S_dd(s_i, opt_num , 5)]; 
         for S_d_n=2:5
              Constraints = [ Constraints , S_aa(s_i, opt_num , S_d_n)<= S_dd(s_i , opt_num , S_d_n)+ S_dd(s_i , opt_num , S_d_n-1)]; 
         end
               Constraints = [ Constraints , sum(S_dd(s_i , opt_num , : ))==1]; 
          
               Constraints = [ Constraints , S_k1(s_i , opt_num ,1)==0.92*S_aa(s_i , opt_num ,3)+0.01*S_bat(s_i , opt_num ,1)];  
               Constraints = [ Constraints , S_bat(s_i , opt_num ,1)==S_w1(s_i , opt_num,1 )+2*S_w1(s_i , opt_num,2 )+4*S_w1(s_i , opt_num,3 )+8*S_w1(s_i, opt_num,4 )];
         for S_w_m=1:4
               Constraints = [ Constraints , S_aa(s_i , opt_num ,3)-(1-S_L1(s_i ,opt_num ,S_w_m ))*M<=S_w1(s_i , opt_num,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 3)];

               Constraints = [ Constraints , 0<=S_w1(s_i , opt_num ,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num ,S_w_m )<= S_L1(s_i ,opt_num ,S_w_m )*M]; 
         end

              Constraints = [ Constraints , S_k1(s_i , opt_num ,2)==0.92*S_aa(s_i , opt_num , 4)+0.01*S_bat(s_i , opt_num ,2)];  
              Constraints = [ Constraints , S_bat(s_i , opt_num ,2)==S_w2(s_i , opt_num,1 )+2*S_w2(s_i , opt_num,2 )+4*S_w2(s_i , opt_num,3 )+8*S_w2(s_i, opt_num,4 )];
         for S_w_m=1:4
              Constraints = [ Constraints , S_aa(s_i , opt_num , 4)-(1-S_L2(s_i ,opt_num ,S_w_m ))*M<=S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 4)]; 
              Constraints = [ Constraints , 0<= S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_L2(s_i ,opt_num ,S_w_m )*M]; 
         end
%              Constraints = [ Constraints ,   1<=S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i ,4 )+16*S_L1(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   1<=S_L2(s_i,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L2(s_i ,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i ,4 )+16*S_L2(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   S_L1(s_i,1 )+2*S_L1(s_i,2 )+4*S_L1(s_i,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=S_L2(s_i,1 )+2*S_L2(s_i,2 )+4*S_L2(s_i,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
   
             Constraints = [ Constraints ,   1<=S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   1<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
      end
     if( s_i == 2 && b1(2 , opt_num)== 1) 
         Constraints = [ Constraints , Q1(s_i , opt_num )== S_v(s_i , opt_num )*Q_sop(s_i , opt_num )];    
         Constraints = [ Constraints , x_ui_square( 25 , opt_num )== 0.9*S_aa(s_i , opt_num ,1)+0.92*S_aa(s_i , opt_num , 2 )+S_k1(s_i , opt_num ,1)+S_k1(s_i , opt_num ,2)+1.08*S_aa(s_i , opt_num , 5 )+1.1*S_aa(s_i , opt_num , 6 )];
         Constraints = [ Constraints , S_v(s_i , opt_num )== S_aa(s_i , opt_num , 1)+ S_aa(s_i, opt_num , 2 )+ S_aa(s_i, opt_num , 3 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i , opt_num )+S_aa(s_i , opt_num , 4 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i , opt_num )...
                                       - S_aa(s_i , opt_num , 5 )- S_aa(s_i , opt_num , 6 )];
                           
         
         Constraints = [ Constraints , 0<=S_aa(s_i , opt_num , : )]; 
         Constraints = [ Constraints , sum(S_aa(s_i , opt_num , : ))==1];                          

         Constraints = [ Constraints , S_aa(s_i , opt_num , 1)<=S_dd(s_i , opt_num , 1)]; 
         Constraints = [ Constraints , S_aa(s_i , opt_num , 6)<=S_dd(s_i, opt_num , 5)]; 
         for S_d_n=2:5
              Constraints = [ Constraints , S_aa(s_i, opt_num , S_d_n)<= S_dd(s_i , opt_num , S_d_n)+ S_dd(s_i , opt_num , S_d_n-1)]; 
         end
               Constraints = [ Constraints , sum(S_dd(s_i , opt_num , : ))==1]; 
          
                Constraints = [ Constraints , S_k1(s_i , opt_num ,1)==0.92*S_aa(s_i , opt_num ,3)+0.01*S_bat(s_i , opt_num ,1)];  
               Constraints = [ Constraints , S_bat(s_i , opt_num ,1)==S_w1(s_i , opt_num,1 )+2*S_w1(s_i , opt_num,2 )+4*S_w1(s_i , opt_num,3 )+8*S_w1(s_i, opt_num,4 )];
         for S_w_m=1:4
               Constraints = [ Constraints , S_aa(s_i , opt_num ,3)-(1-S_L1(s_i ,opt_num ,S_w_m ))*M<=S_w1(s_i , opt_num,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 3)];

               Constraints = [ Constraints , 0<=S_w1(s_i , opt_num ,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num ,S_w_m )<= S_L1(s_i ,opt_num ,S_w_m )*M]; 
         end

              Constraints = [ Constraints , S_k1(s_i , opt_num ,2)==0.92*S_aa(s_i , opt_num , 4)+0.01*S_bat(s_i , opt_num ,2)];  
              Constraints = [ Constraints , S_bat(s_i , opt_num ,2)==S_w2(s_i , opt_num,1 )+2*S_w2(s_i , opt_num,2 )+4*S_w2(s_i , opt_num,3 )+8*S_w2(s_i, opt_num,4 )];
         for S_w_m=1:4
              Constraints = [ Constraints , S_aa(s_i , opt_num , 4)-(1-S_L2(s_i ,opt_num ,S_w_m ))*M<=S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 4)]; 
              Constraints = [ Constraints , 0<= S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_L2(s_i ,opt_num ,S_w_m )*M]; 
         end
%              Constraints = [ Constraints ,   1<=S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i ,4 )+16*S_L1(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   1<=S_L2(s_i,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L2(s_i ,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i ,4 )+16*S_L2(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   S_L1(s_i,1 )+2*S_L1(s_i,2 )+4*S_L1(s_i,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=S_L2(s_i,1 )+2*S_L2(s_i,2 )+4*S_L2(s_i,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
   
             Constraints = [ Constraints ,   1<=S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   1<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
      end
     if( s_i == 3 && b2(1 , opt_num)== 1 )  
         Constraints = [ Constraints , Q2(s_i-2  , opt_num )== S_v(s_i , opt_num )*Q_sop(s_i , opt_num )];    
         Constraints = [ Constraints , x_ui_square( 22 , opt_num )== 0.9*S_aa(s_i , opt_num ,1)+0.92*S_aa(s_i , opt_num , 2 )+S_k1(s_i , opt_num ,1)+S_k1(s_i , opt_num ,2)+1.08*S_aa(s_i , opt_num , 5 )+1.1*S_aa(s_i , opt_num , 6 )];
         Constraints = [ Constraints , S_v(s_i , opt_num )== S_aa(s_i , opt_num , 1)+ S_aa(s_i, opt_num , 2 )+ S_aa(s_i, opt_num , 3 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i , opt_num )+S_aa(s_i , opt_num , 4 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i ,opt_num)...
                                       - S_aa(s_i , opt_num , 5 )- S_aa(s_i , opt_num , 6 )];
                                   
         Constraints = [ Constraints , 0<=S_aa(s_i , opt_num , : )]; 
         Constraints = [ Constraints , sum(S_aa(s_i , opt_num , : ))==1];                          

         Constraints = [ Constraints , S_aa(s_i , opt_num , 1)<=S_dd(s_i , opt_num , 1)]; 
         Constraints = [ Constraints , S_aa(s_i , opt_num ,6)<=S_dd(s_i, opt_num , 5)]; 
         for S_d_n=2:5
              Constraints = [ Constraints , S_aa(s_i, opt_num , S_d_n)<= S_dd(s_i , opt_num , S_d_n)+ S_dd(s_i , opt_num , S_d_n-1)]; 
         end
               Constraints = [ Constraints , sum(S_dd(s_i , opt_num , : ))==1]; 
          
                 Constraints = [ Constraints , S_k1(s_i , opt_num ,1)==0.92*S_aa(s_i , opt_num ,3)+0.01*S_bat(s_i , opt_num ,1)];  
               Constraints = [ Constraints , S_bat(s_i , opt_num ,1)==S_w1(s_i , opt_num,1 )+2*S_w1(s_i , opt_num,2 )+4*S_w1(s_i , opt_num,3 )+8*S_w1(s_i, opt_num,4 )];
         for S_w_m=1:4
               Constraints = [ Constraints , S_aa(s_i , opt_num ,3)-(1-S_L1(s_i ,opt_num ,S_w_m ))*M<=S_w1(s_i , opt_num,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 3)];

               Constraints = [ Constraints , 0<=S_w1(s_i , opt_num ,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num ,S_w_m )<= S_L1(s_i ,opt_num ,S_w_m )*M]; 
         end

              Constraints = [ Constraints , S_k1(s_i , opt_num ,2)==0.92*S_aa(s_i , opt_num , 4)+0.01*S_bat(s_i , opt_num ,2)];  
              Constraints = [ Constraints , S_bat(s_i , opt_num ,2)==S_w2(s_i , opt_num,1 )+2*S_w2(s_i , opt_num,2 )+4*S_w2(s_i , opt_num,3 )+8*S_w2(s_i, opt_num,4 )];
         for S_w_m=1:4
              Constraints = [ Constraints , S_aa(s_i , opt_num , 4)-(1-S_L2(s_i ,opt_num ,S_w_m ))*M<=S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 4)]; 
              Constraints = [ Constraints , 0<= S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_L2(s_i ,opt_num ,S_w_m )*M]; 
         end
%              Constraints = [ Constraints ,   1<=S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i ,4 )+16*S_L1(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   1<=S_L2(s_i,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L2(s_i ,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i ,4 )+16*S_L2(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   S_L1(s_i,1 )+2*S_L1(s_i,2 )+4*S_L1(s_i,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=S_L2(s_i,1 )+2*S_L2(s_i,2 )+4*S_L2(s_i,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
   
             Constraints = [ Constraints ,   1<=S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   1<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
      end
     if( s_i == 4 && b2(2 , opt_num)== 1)  
         Constraints = [ Constraints , Q2(s_i-2 , opt_num )== S_v(s_i , opt_num )*Q_sop(s_i , opt_num )];    
         Constraints = [ Constraints , x_ui_square( 29 , opt_num )== 0.9*S_aa(s_i , opt_num ,1)+0.92*S_aa(s_i , opt_num , 2 )+S_k1(s_i , opt_num ,1)+S_k1(s_i , opt_num ,2)+1.08*S_aa(s_i , opt_num , 5 )+1.1*S_aa(s_i , opt_num , 6 )];
         Constraints = [ Constraints , S_v(s_i , opt_num )== S_aa(s_i , opt_num , 1)+ S_aa(s_i, opt_num , 2 )+ S_aa(s_i, opt_num , 3 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i , opt_num )+S_aa(s_i , opt_num , 4 )*Qsop_ref(s_i , opt_num )/Q_sop(s_i , opt_num)...
                                       - S_aa(s_i , opt_num , 5 )- S_aa(s_i , opt_num , 6 )];
                                   
         Constraints = [ Constraints , 0<=S_aa(s_i , opt_num , : )]; 
         Constraints = [ Constraints , sum(S_aa(s_i , opt_num , : ))==1];                          

         Constraints = [ Constraints , S_aa(s_i , opt_num , 1)<=S_dd(s_i , opt_num , 1)]; 
         Constraints = [ Constraints , S_aa(s_i , opt_num , 6)<=S_dd(s_i, opt_num , 5)]; 
         for S_d_n=2:5
              Constraints = [ Constraints , S_aa(s_i, opt_num , S_d_n)<= S_dd(s_i , opt_num , S_d_n)+ S_dd(s_i , opt_num , S_d_n-1)]; 
         end
               Constraints = [ Constraints , sum(S_dd(s_i , opt_num , : ))==1]; 
               
               Constraints = [ Constraints , S_k1(s_i , opt_num ,1)==0.92*S_aa(s_i , opt_num ,3)+0.01*S_bat(s_i , opt_num ,1)];  
               Constraints = [ Constraints , S_bat(s_i , opt_num ,1)==S_w1(s_i , opt_num,1 )+2*S_w1(s_i , opt_num,2 )+4*S_w1(s_i , opt_num,3 )+8*S_w1(s_i, opt_num,4 )];
         for S_w_m=1:4
               Constraints = [ Constraints , S_aa(s_i , opt_num ,3)-(1-S_L1(s_i ,opt_num ,S_w_m ))*M<=S_w1(s_i , opt_num,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 3)];

               Constraints = [ Constraints , 0<=S_w1(s_i , opt_num ,S_w_m )]; 
               Constraints = [ Constraints , S_w1(s_i , opt_num ,S_w_m )<= S_L1(s_i ,opt_num ,S_w_m )*M]; 
         end

              Constraints = [ Constraints , S_k1(s_i , opt_num ,2)==0.92*S_aa(s_i , opt_num , 4)+0.01*S_bat(s_i , opt_num ,2)];  
              Constraints = [ Constraints , S_bat(s_i , opt_num ,2)==S_w2(s_i , opt_num,1 )+2*S_w2(s_i , opt_num,2 )+4*S_w2(s_i , opt_num,3 )+8*S_w2(s_i, opt_num,4 )];
         for S_w_m=1:4
              Constraints = [ Constraints , S_aa(s_i , opt_num , 4)-(1-S_L2(s_i ,opt_num ,S_w_m ))*M<=S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_aa(s_i , opt_num , 4)]; 
              Constraints = [ Constraints , 0<= S_w2(s_i , opt_num,S_w_m )]; 
              Constraints = [ Constraints , S_w2(s_i , opt_num,S_w_m )<= S_L2(s_i ,opt_num ,S_w_m )*M]; 
         end
%              Constraints = [ Constraints ,   1<=S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i ,4 )+16*S_L1(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L1(s_i ,1 )+2*S_L1(s_i ,2 )+4*S_L1(s_i ,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   1<=S_L2(s_i,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
%              Constraints = [ Constraints ,   S_L2(s_i ,1 )+2*S_L2(s_i ,2 )+4*S_L2(s_i ,3 )+8*S_L2(s_i ,4 )+16*S_L2(s_i ,5 )<=19]; 
%              Constraints = [ Constraints ,   S_L1(s_i,1 )+2*S_L1(s_i,2 )+4*S_L1(s_i,3 )+8*S_L1(s_i,4 )+16*S_L1(s_i ,5 )<=S_L2(s_i,1 )+2*S_L2(s_i,2 )+4*S_L2(s_i,3 )+8*S_L2(s_i,4 )+16*S_L2(s_i ,5 )]; 
   
             Constraints = [ Constraints ,   1<=S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   1<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
             Constraints = [ Constraints ,   S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )<=15]; 
             Constraints = [ Constraints ,   S_L1(s_i ,opt_num ,1 )+2*S_L1(s_i ,opt_num ,2 )+4*S_L1(s_i ,opt_num ,3 )+8*S_L1(s_i ,opt_num ,4 )<=S_L2(s_i,opt_num ,1 )+2*S_L2(s_i ,opt_num ,2 )+4*S_L2(s_i ,opt_num ,3 )+8*S_L2(s_i,opt_num ,4 )]; 
      end
end

%% 光伏机组功率削减约束
% Constraints = [ Constraints , 0<= PV_cut( 1 , opt_num )<= p_Solar(pv1 , opt_num )];
% Constraints = [ Constraints , 0<= PV_cut( 2 , opt_num )<= p_Solar(pv2 , opt_num )];
% Constraints = [ Constraints , 0<= PV_cut( 3 , opt_num )<= p_Solar(pv3 , opt_num )];

Constraints = [ Constraints , P_pv( 1 , opt_num )== p_Solar(pv1 , opt_num )-PV_cut( 1 , opt_num )];
Constraints = [ Constraints , P_pv( 2 , opt_num )== p_Solar(pv2 , opt_num )-PV_cut( 2 , opt_num )];
Constraints = [ Constraints , P_pv( 3 , opt_num )== p_Solar(pv3 , opt_num )-PV_cut( 3 , opt_num )];
Constraints = [ Constraints , P_pv( 4 , opt_num )== p_Solar(pv4 , opt_num )-PV_cut( 4 , opt_num )];
Constraints = [ Constraints , P_pv( 5 , opt_num )== p_Solar(pv5 , opt_num )-PV_cut( 5 , opt_num )];
Constraints = [ Constraints , P_pv( 6 , opt_num ) == p_Solar(pv6 , opt_num )-PV_cut( 6 , opt_num )];
%% 变压器1-2节点电压变比约束
Constraints = [ Constraints , x_ui_square( 2 , opt_num )  == ( k_ij0 + Ktij( opt_num ) * delta_kij ) * x_ui_square( 1 , opt_num )];

Constraints = [ Constraints , Aux( : , opt_num ) >= x_ui_square( : , opt_num )  - 1.03 ];
Constraints = [ Constraints , Aux( : , opt_num ) >= -x_ui_square( : , opt_num )  + 0.97 ];
Constraints = [ Constraints , Aux( : , opt_num ) >= 0 ];
toc
end
%% 初始运行点ES-SOP运行约束
for opt=1:N
for rr=1:2
%ACDC损耗
Constraints = [ Constraints , P1_L( rr , opt) == C_a0*b1( rr , opt ) + C_a1*k11( rr , opt) + C_a2*k12(rr , opt) ];
Constraints = [ Constraints , P2_L( rr , opt ) == C_a0*b2( rr , opt ) + C_a1*k21( rr, opt ) + C_a2*k22( rr , opt) ];
%ACDC辅助变量定义松弛
Constraints = [ Constraints , norm( [P1( rr , opt) Q1( rr , opt )])<= k11( rr , opt) ];
Constraints = [ Constraints , norm( [2*P1( rr , opt ) 2*Q1( rr , opt) (1-k12( rr , opt))])<= 1+ k12( rr , opt)];
Constraints = [ Constraints , norm( [P2( rr , opt) Q2( rr , opt)])<= k21( rr , opt) ];
Constraints = [ Constraints , norm( [2*P2( rr , opt) 2*Q2( rr , opt ) (1-k22( rr , opt))])<= 1+ k22( rr , opt )];
%ACDC切断控制
Constraints = [ Constraints ,k11( rr , opt )<=b1( rr , opt)*S_acdc ];
Constraints = [ Constraints ,k12( rr , opt )<=b1( rr , opt )*(S_acdc)^2 ];
Constraints = [ Constraints ,k21( rr , opt)<=b2( rr , opt )*S_acdc ];
Constraints = [ Constraints ,k22( rr , opt )<=b2( rr , opt )*(S_acdc)^2 ];

%容量约束
Constraints = [ Constraints , norm( [P1( rr , opt) Q1( rr , opt )])<= S_acdc ];
Constraints = [ Constraints , norm( [P2( rr , opt) Q2( rr , opt )])<= S_acdc ];
end
end
%----------------目标函数求和------------------------------------------
f1 = sum(f);

%----------------求解器及求解设置------------------------------------------
% options = sdpsettings('verbose',1,'solver','gurobi');%gurobi
options = sdpsettings('verbose',1,'solver','gurobi' , 'gurobi.MIPGap' , 0.005);%gurobi

sol = solvesdp(Constraints,f1,options);
% --------------------------Analyze error flags----------------------------
if (sol.problem == 0)
 % Extract and display value
    result = double(f1)
else
    
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
%% 数据保存
%------------------------电压、电流、支路有功无功--------------------------
x_ui_squaret(:,t_1:t_2)=x_ui_square(:,1:N);
x_pijt(:,t_1:t_2)=x_pij(:,1:N); 
x_qijt(:,t_1:t_2)=x_qij(:,1:N);
%------------------------下垂控制优化参数--------------------------
% for q_i=1:q_n
% v_min(q_i,t)= 0.92+0.01*(L1(q_i,1 )+2*L1(q_i,2 )+4*L1(q_i,3 )+8*L1(q_i,4 )+16*L1(q_i,5 ));
% v_max(q_i,t)=0.92+0.01*(L2(q_i,1 )+2*L2(q_i,2 )+4*L2(q_i,3 )+8*L2(q_i,4 )+16*L2(q_i,5 )); 
% end
% for q_i=1:4
% S_v_min(q_i,t)= 0.92+0.01*(S_L1(q_i,1 )+2*S_L1(q_i,2 )+4*S_L1(q_i,3 )+8*S_L1(q_i,4 )+16*S_L1(q_i,5 ));
% S_v_max(q_i,t)=0.92+0.01*(S_L2(q_i,1 )+2*S_L2(q_i,2 )+4*S_L2(q_i,3 )+8*S_L2(q_i,4 )+16*S_L2(q_i,5 )); 
% end
for q_i=1:q_n
v_min(q_i,t)= 0.92+0.01*(L1(q_i,1,1 )+2*L1(q_i,1,2 )+4*L1(q_i,1,3 )+8*L1(q_i,1,4 ));
v_max(q_i,t)=0.92+0.01*(L2(q_i,1,1 )+2*L2(q_i,1,2 )+4*L2(q_i,1,3 )+8*L2(q_i,1,4 )); 
end
for q_i=1:4
S_v_min(q_i,t)= 0.92+0.01*(S_L1(q_i,1,1 )+2*S_L1(q_i,1,2 )+4*S_L1(q_i,1,3 )+8*S_L1(q_i,1,4 ));
S_v_max(q_i,t)=0.92+0.01*(S_L2(q_i,1,1 )+2*S_L2(q_i,1,2 )+4*S_L2(q_i,1,3 )+8*S_L2(q_i,1,4 )); 
end
%------------------------下垂控制无功标幺值--------------------------
v_f(:,t_1:t_2)=v(:,1:N); 
S_v_f(:,t_1:t_2)=S_v(:,1:N);
%-----------------------损耗-----------------------------------
f_startt(:,t_1:t_2)=f_start(:,1:N);
f_essopt(:,t_1:t_2)=f_essop(:,1:N);
f_PV_cutt(:,t_1:t_2)=f_PV_cut(:,1:N); 
%-----------------------ES-SOP-----------------------------------
%AC/DC无功功率
Q1t(:,t_1:t_2)=Q1(:,1:N); 
Q2t(:,t_1:t_2)=Q2(:,1:N); 
%-----------------------光伏无功输出-----------------------------------
q_Solart(:,t_1:t_2)=q_Solar(:,1:N);
PV_cutt(:,t_1:t_2)=PV_cut(:,1:N);

end
disp(['求解完成',num2str(etime(clock,t1))]);

%% 数据读取
%------------------------下垂控制优化参数--------------------------
v_min_data20 = double(v_min);
v_max_data20= double(v_max);
v_f_data20=double(v_f);

S_v_min_data20 = double(S_v_min);
S_v_max_data20 = double(S_v_max);
S_v_f_data20=double(S_v_f);
%------------------------电压、支路有功无功--------------------------
x_ui_square_data20 = double((x_ui_squaret));
x_pij_data20= double(x_pijt);
x_qij_data20= double(x_qijt);

%-----------------------损耗-----------------------------------
f_start_data20=double(f_startt);
f_essop_data20=double(f_essopt);
f_PV_cut_data20=double(f_PV_cutt);
%-----------------------ES-SOP-----------------------------------
%AC/DC无功功率
Q1_data20=double(Q1t);
Q2_data20=double(Q2t);
%-----------------------电压偏差-----------------------------------
aux20=sum(sum(abs(x_ui_square_data20-1)))*1/4;
%-----------------------光伏无功输出-----------------------------------
q_Solar_data20=double(q_Solart);
PV_cut_data20=double(PV_cutt);
%% 存储数据
save stage2_2_WithoutVSM_随机性 x_ui_square_data20  x_pij_data20 x_qij_data20...
f_start_data20  f_PV_cut_data20 f_essop_data20...
Q1_data20 Q2_data20 ...
aux20...
q_Solar_data20 PV_cut_data20...
v_min_data20 v_max_data20 S_v_min_data20 S_v_max_data20 v_f_data20 S_v_f_data20;

save 下垂参数_WithoutVSM v_min_data20 v_max_data20 S_v_min_data20 S_v_max_data20;

% a1=sum(f_start_data)
% a2=sum(f_essop_data)
 
% u1_min=min(x_ui_square_data');
% u1_max=max(x_ui_square_data');

% u2_min=min(x_ui_square_data');
% u2_max=max(x_ui_square_data');
 
% u3_min=min(x_ui_square_data');
% u3_max=max(x_ui_square_data');

% save 每个节点最大最小电压标幺值_有vsm u1_min u1_max u2_min u2_max u3_min u3_max;
% save dianya u1_min u1_max u2_min u2_max;

