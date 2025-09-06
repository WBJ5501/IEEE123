clear all
clc

run ieee_33_node_system.m;
run DG_Load.m;
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
result = zeros(1,N);
%------------------电压上下限--------------------------------------
a1=0.95;
a2=1.05;
%----------------变压器抽头及电容器--------------------------------
delta_T = 1;%hour
N =24;%时间
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
% Q_sop = 2*0.6;

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

%----------------负荷光伏数据及恒定功率因数--------------------------------
p_Solar = zeros(33,N);
p_Wind  = zeros(33,N);
p_Load  = zeros(33,N);
q_Load  = zeros(33,N);

tan_theta = q_load ./ p_load;%恒定功率因数

for a = 1 : N
    p_Solar( : , a ) = Solar_radio * Solar_data( a );
    p_Load ( : , a ) = Load_radio  *  Load_data( a );
    q_Load ( : , a ) = ( p_Load ( : , a ) ./ p_load ) .* ( q_load )  ;
end


%% 定义初始运行点决策变量
q_Solar =sdpvar(q_n, N, 'full');
%-----------------光伏削减量------------------------------------
PV_cut=sdpvar(q_n, N, 'full');
P_pv=sdpvar(q_n, N, 'full');

%----------------电压、支路有功、支路无功--------------------------------
x_ui_square = sdpvar(33, N, 'full');
x_pij = sdpvar(32, N, 'full');
x_qij = sdpvar(32, N, 'full');

%----------------变压器抽头--------------------------------
btij = binvar(K_max*2+1, N, 'full');

Ktij = intvar(1, N+1, 'full');
Ktij1 = intvar(1, N, 'full');
Ktij2 = intvar(1, N, 'full');

%----------------电容器--------------------------------
qc=3;
NtCB = intvar(qc, N+1, 'full');
NtCB1 = intvar(qc, N, 'full');
NtCB2 = intvar(qc, N, 'full');

%--------------------1组ESSOP/双馈线-------------------------------------
sv=2;
%VSC接通/断开变量
b1 = binvar(sv, N, 'full');
b2 = binvar(sv, N, 'full');
b3 = binvar(sv-1, N, 'full');
%AC/DC视在功率辅助变量/非负
k11= sdpvar(sv, N, 'full');
k12= sdpvar(sv, N, 'full');
k21= sdpvar(sv, N, 'full');
k22= sdpvar(sv, N, 'full');
%四个子系统损耗
P1_L=sdpvar(sv, N, 'full');
P2_L=sdpvar(sv, N, 'full');
P3_L=sdpvar(sv-1, N, 'full');
P4_L=sdpvar(sv-1, N, 'full');
%四个子系统有功功率
P1=sdpvar(sv, N, 'full');
P2=sdpvar(sv, N, 'full');
P3=sdpvar(sv-1, N, 'full');
P4=sdpvar(sv-1, N, 'full');
%AC/DC无功功率
Q1=sdpvar(sv, N, 'full');
Q2=sdpvar(sv, N, 'full');
%DC/DC损耗中的绝对值项和平方项
P3_abs=sdpvar(sv-1, N, 'full');
k3=sdpvar(sv-1, N, 'full');
%ESS损耗中的平方项和荷电状态
k4=sdpvar(sv-1, N, 'full');    
Soc_ess=sdpvar(sv-1, N+1, 'full');

%----------------损耗-----------------------------
P_loss=sdpvar(32, N, 'full'); 

%% 定义约束条件
Constraints = [];
% Constraints = [ Constraints , q_Solar== 0 ];
% Constraints = [ Constraints , b1 == 0 ];
% Constraints = [ Constraints , b2 == 0 ];
% Constraints = [ Constraints , b3 == 0 ];

for opt_num = 1: (N)
     
    tic
%% 目标函数

        f(opt_num)=1000*0.08*(sum(P_loss( : , opt_num ))+sum( P1_L( : , opt_num )+ P2_L( : , opt_num ) )+ P3_L( : , opt_num )+ P4_L( : , opt_num ))...
                  + 1.4 * ( Ktij1(opt_num) + Ktij2(opt_num) ) + 0.24 * sum( NtCB1(:,opt_num) + NtCB2(:,opt_num))+ 1000*0.64*sum(PV_cut( : , opt_num ));   


        f_start(opt_num)= sum(P_loss( : , opt_num ));
        f_essop(opt_num)= sum( P1_L( : , opt_num )+ P2_L( : , opt_num ))+ P3_L( : , opt_num )+ P4_L( : , opt_num ) ;
        f_oltc(opt_num)= sum( Ktij1(opt_num) + Ktij2(opt_num) );
        f_ntcb(opt_num)= sum( NtCB1(:,opt_num) + NtCB2(:,opt_num));
        f_PV_cut(opt_num)= sum(PV_cut( : , opt_num ));

        
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
                                    + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num )== 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                    - q_Load( k , opt_num )== 0 ];                        
  
    elseif( k == 33 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num ) )...
                                    + p_Solar( k , opt_num ) + p_Wind( k , opt_num )- p_Load( k , opt_num ) == 0 ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                    - q_Load( k , opt_num  )== 0 ];
                                                                                             
    elseif( k == 6 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )+ NtCB(1,opt_num + 1  )*qCB   == sum( x_qij( node_out , opt_num ) ) ];
                                  
    elseif( k == 13 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num ) )...
                                      - q_Load(k , opt_num )+ NtCB( 2,opt_num + 1  )*qCB   == sum( x_qij( node_out , opt_num ) ) ];                                                             
                                  
    elseif( k == 31 )
        Constraints = [ Constraints , sum( x_pij( node_in , opt_num )  )...
                                      + p_Solar( k , opt_num ) + p_Wind( k , opt_num ) - p_Load( k , opt_num ) == sum( x_pij( node_out , opt_num ) ) ];
        Constraints = [ Constraints , sum( x_qij( node_in , opt_num )  )...
                                      - q_Load(k , opt_num )  + NtCB( 3,opt_num + 1  )*qCB  == sum( x_qij( node_out , opt_num ) ) ];
                                                      
                                         
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

%% 初始运行点ES-SOP运行约束
for rr=1:2
%ACDC损耗
Constraints = [ Constraints , P1_L( rr , opt_num ) == C_a0*b1( rr , opt_num ) + C_a1*k11( rr , opt_num ) + C_a2*k12(rr , opt_num ) ];
Constraints = [ Constraints , P2_L( rr , opt_num ) == C_a0*b2( rr , opt_num ) + C_a1*k21( rr, opt_num ) + C_a2*k22( rr , opt_num ) ];
%ACDC辅助变量定义松弛
Constraints = [ Constraints , norm( [P1( rr , opt_num ) Q1( rr , opt_num )])<= k11( rr , opt_num ) ];
Constraints = [ Constraints , norm( [2*P1( rr , opt_num ) 2*Q1( rr , opt_num ) (1-k12( rr , opt_num ))])<= 1+ k12( rr , opt_num )];
Constraints = [ Constraints , norm( [P2( rr , opt_num ) Q2( rr , opt_num )])<= k21( rr , opt_num ) ];
Constraints = [ Constraints , norm( [2*P2( rr , opt_num ) 2*Q2( rr , opt_num ) (1-k22( rr , opt_num ))])<= 1+ k22( rr , opt_num )];
%ACDC切断控制
Constraints = [ Constraints ,k11( rr , opt_num )<=b1( rr , opt_num )*S_acdc ];
Constraints = [ Constraints ,k12( rr , opt_num )<=b1( rr , opt_num )*(S_acdc)^2 ];
Constraints = [ Constraints ,k21( rr , opt_num )<=b2( rr , opt_num )*S_acdc ];
Constraints = [ Constraints ,k22( rr , opt_num )<=b2( rr , opt_num )*(S_acdc)^2 ];

%容量约束
Constraints = [ Constraints , norm( [P1( rr , opt_num ) Q1( rr , opt_num )])<= S_acdc ];
Constraints = [ Constraints , norm( [P2( rr , opt_num ) Q2( rr , opt_num )])<= S_acdc ];
end


%DCDC损耗
Constraints = [ Constraints , P3_L( 1 , opt_num ) == C_d0*b3( 1 , opt_num ) + C_d1*P3_abs( 1 , opt_num ) + C_d2*k3(1 , opt_num ) ];
%DCDC辅助变量定义松弛
Constraints = [ Constraints ,P3( 1 , opt_num )<=P3_abs( 1 , opt_num ) ];
Constraints = [ Constraints ,-P3( 1 , opt_num )<=P3_abs( 1 , opt_num ) ];
Constraints = [ Constraints , norm( [2*P3( 1 , opt_num ) (1-k3( 1 , opt_num ))])<= 1+ k3( 1 , opt_num )];
%DCDC切断控制
Constraints = [ Constraints ,P3_abs( 1 , opt_num )<=b3( 1 , opt_num )*S_dcdc ];
Constraints = [ Constraints ,k3( 1 , opt_num )<=b3( 1 , opt_num )*(S_dcdc)^2 ];

%ESS损耗
Constraints = [ Constraints , P4_L( 1 , opt_num ) == (C_aux + C_ES*k4( 1 , opt_num ))];
Constraints = [ Constraints , norm( [2*P4( 1 , opt_num ) (1-k4( 1 , opt_num ))])<= 1+ k4( 1 , opt_num )];


%储能运行约束
P_max=0.6;
S_min=0.2;
S_max=1;
Soc_ess( 1  , 1 )=0.5;
Constraints = [ Constraints , Soc_ess( 1 , opt_num+1 ) == Soc_ess( 1 , opt_num )-( P4_L( 1 , opt_num ) + P4( 1 , opt_num ) )];
Constraints = [ Constraints ,-P_max <= P4( 1 , opt_num ) ];
Constraints = [ Constraints , P4( 1 , opt_num ) <= P_max ];
Constraints = [ Constraints ,S_min <= Soc_ess( 1 , opt_num )  ];
Constraints = [ Constraints ,Soc_ess( 1 , opt_num ) <= S_max ];
Constraints = [ Constraints ,Soc_ess( 1  , 1 ) == Soc_ess( 1 , 25 )];

%容量约束
Constraints = [ Constraints ,-S_dcdc <= P3( 1 , opt_num ) ];
Constraints = [ Constraints , P3( 1 , opt_num ) <=S_dcdc ];

%有功平衡
Constraints = [ Constraints , P3( 1 , opt_num ) == -(P4( 1 , opt_num )-P4_L( 1 , opt_num ))];

Constraints = [ Constraints , sum(P1( : , opt_num ) + P2( : , opt_num )+ P1_L( : , opt_num )+ P2_L( : , opt_num ))+ P3( 1 , opt_num ) + P3_L( 1 , opt_num )==0];

%多馈线约束
Constraints = [ Constraints , sum(b1( : , opt_num )) <= 1 ];
Constraints = [ Constraints , sum(b1( 1 , opt_num )+ b2( 2 , opt_num )) <= 1 ];
Constraints = [ Constraints , sum(b2( : , opt_num ))  <= 1 ];
Constraints = [ Constraints , sum(b1( 2 , opt_num )+ b2( 1 , opt_num ))<= 1 ];

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

%% 光伏机组功率削减约束
Constraints = [ Constraints , 0<= PV_cut( 1 , opt_num )<= p_Solar(pv1 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 2 , opt_num )<= p_Solar(pv2 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 3 , opt_num )<= p_Solar(pv3 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 4 , opt_num )<= p_Solar(pv4 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 5 , opt_num )<= p_Solar(pv5 , opt_num )];
Constraints = [ Constraints , 0<= PV_cut( 6 , opt_num )<= p_Solar(pv6 , opt_num )];


Constraints = [ Constraints , P_pv( 1 , opt_num )== p_Solar(pv1 , opt_num )-PV_cut( 1 , opt_num )];
Constraints = [ Constraints , P_pv( 2 , opt_num )== p_Solar(pv2 , opt_num )-PV_cut( 2 , opt_num )];
Constraints = [ Constraints , P_pv( 3 , opt_num )== p_Solar(pv3 , opt_num )-PV_cut( 3 , opt_num )];
Constraints = [ Constraints , P_pv( 4 , opt_num )== p_Solar(pv4 , opt_num )-PV_cut( 4 , opt_num )];
Constraints = [ Constraints , P_pv( 5 , opt_num )== p_Solar(pv5 , opt_num )-PV_cut( 5 , opt_num )];
Constraints = [ Constraints , P_pv( 6 , opt_num ) == p_Solar(pv6 , opt_num )-PV_cut( 6 , opt_num )];
%% 变压器1-2节点电压变比约束
Constraints = [ Constraints , Ktij( opt_num ) == (kk - K_max) * btij( : , opt_num ) ];
Constraints = [ Constraints , sum( btij( : , opt_num ) ) == 1 ];
Constraints = [ Constraints , x_ui_square( 2 , opt_num )  == ( k_ij0 + Ktij( opt_num ) * delta_kij ) * x_ui_square( 1 , opt_num )];

toc
end

%% OLTC约束
Constraints = [ Constraints , Ktij( 1 ) == 0 ];
Constraints = [ Constraints , sum( Ktij1 + Ktij2 ) <= (K_max-1) ];
for k = 2 : N+1
        Constraints = [ Constraints , Ktij( k ) - Ktij( k - 1 ) ==  Ktij1(k-1) - Ktij2(k-1)];
end
Constraints = [ Constraints , Ktij1 >= 0 ];
Constraints = [ Constraints , Ktij2 >= 0 ];
Constraints = [ Constraints , -K_max <= Ktij <= K_max ];

%% CB约束
for n=1:qc
    
Constraints = [ Constraints , NtCB(n, 1 ) == 1 ];
Constraints = [ Constraints , sum( NtCB1(n, : ) + NtCB2(n, : ) ) <= 4 ];

for k = 2 : N+1
        Constraints = [ Constraints , NtCB(n, k ) - NtCB(n, k - 1 ) == NtCB1(n,k-1) - NtCB2(n,k-1) ];
end

Constraints = [ Constraints , NtCB1(n, : ) >= 0 ];
Constraints = [ Constraints , NtCB2(n, : ) >= 0 ];
Constraints = [ Constraints , NtCB(n, : ) <= 5 ];
end

%----------------目标函数求和------------------------------------------
f1 = sum(f);

%----------------求解器及求解设置------------------------------------------
options = sdpsettings('verbose',1,'solver','gurobi');%gurobi
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
disp(['求解完成',num2str(etime(clock,t1))]);



%% 数据读取-VSM
%------------------------电压、支路有功无功--------------------------
x_ui_square_data = double((x_ui_square));
x_pij_data = double(x_pij);
x_qij_data = double(x_qij);
%-----------------------临界节点电压-----------------------------------

f_data=double(f);%目标函数
f_start_data=double(f_start);%网损
f_essop_data=double(f_essop);%essop
f_oltc_data=double(f_oltc);
f_ntcb_data=double(f_ntcb);
f_PV_cut_data=double(f_PV_cut);


%-----------------------ES-SOP-----------------------------------
 %VSC接通/断开变量
b1_data = double(b1);
b2_data = double(b2);
b3_data = double(b3);

%AC/DC视在功率辅助变量/非负
k11_data= double(k11);
k12_data= double(k12);
k21_data= double(k21);
k22_data= double(k22);
%四个子系统损耗
P1_L_data=double(P1_L);
P2_L_data=double(P2_L);
P3_L_data=double(P3_L);
P4_L_data=double(P4_L);
%四个子系统有功功率
P1_data=double(P1);
P2_data=double(P2);
P3_data=double(P3);
P4_data=double(P4);
%AC/DC无功功率
Q1_data=double(Q1);
Q2_data=double(Q2);
%DC/DC损耗中的绝对值项和平方项
P3_abs_data=double(P3_abs);
k3_data=double(k3);
%ESS损耗中的平方项和荷电状态
k4_data=double(k4);  
Soc_ess_data=double(Soc_ess);
%-----------------------OLTC-----------------------------------
Ktij_data = double(Ktij);
%-----------------------CBs-----------------------------------
NtCB_data = double(NtCB);

%-----------------------电压偏差-----------------------------------
aux=sum(sum(abs(x_ui_square_data-1)));

%-----------------------光伏无功输出-----------------------------------
q_Solar_data=double(q_Solar);
PV_cut_data=double(PV_cut);









%% 存储数据
save stage1_duokuixian_withoutvsm x_ui_square_data  x_pij_data x_qij_data...
f_start_data f_data f_essop_data f_oltc_data  f_ntcb_data f_PV_cut_data...
b1_data b2_data  b3_data ...  
k11_data k21_data...
P1_data P2_data P3_data P4_data...
P1_L_data P2_L_data P3_L_data P4_L_data...
Q1_data Q2_data ...
Ktij_data...
NtCB_data ...
aux ...
q_Solar_data PV_cut_data;
    
 
% a1=sum(f_start_data)
% a2=sum(f_essop_data)
% a3=sum(L_f_start_data)
% a4=sum(L_f_essop_data)




