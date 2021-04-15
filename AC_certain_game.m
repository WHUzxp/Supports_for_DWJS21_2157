%确定性优化，按照经验分布，交流潮流，拟合电价
%第九个时间断面开始，08:00-次日08:00
clear
clc
%% 基础数据
N=50;%电动汽车数量
%电价数据(08:00-09:00为时段1)
load DLMP_data
load DLMP_data_RT
a_data=[a_data(9:24);a_data(1:8)];b_data=[b_data(9:24);b_data(1:8)];price_RT=[price_RT(:,9:24),price_RT(:,1:8)];
%电动汽车聚类数据(08:00-09:00为时段1)
EVdata=[6,40,15,10,24,105;6,32,16,2,9,132;3,24,8,11,23,89;6,24,12,13,22,97;6,40,25,1,8,101;6,40,16,12,23,175;3,24,8,10,24,35;10,40,20,2,8,88;10,40,18,11,24,112;10,64,25,11,23,66];
%不同类型电动汽车的比例
ratio=EVdata(:,6)/1000;
%电动汽车状态矩阵
X=zeros(24,10);Y=zeros(24,10);
for i=1:10
    X(EVdata(i,4):EVdata(i,5),i)=1;%电动汽车停泊矩阵
end
%% 建模
price_EV=sdpvar(24,1);
pch=sdpvar(24,10);%电动汽车充电
pdis=sdpvar(24,10);%电动汽车放电
S_EV=sdpvar(24,10);%电动汽车电量状态
obj_EV=N*price_EV'*(pch-pdis)*ratio;%电动汽车目标函数
C_EV=[0<=pch<=ones(24,1)*EVdata(:,1)'.*X,0<=pdis<=ones(24,1)*EVdata(:,1)'.*X,
    0.2*ones(24,1)*EVdata(:,2)'.*X<=S_EV<=0.95*ones(24,1)*EVdata(:,2)'.*X];%电动汽车边界约束条件
for n=1:10
    C_EV=[C_EV,S_EV(EVdata(n,4)+1:EVdata(n,5),n)==S_EV(EVdata(n,4):EVdata(n,5)-1,n)+0.95*pch(EVdata(n,4)+1:EVdata(n,5),n)-pdis(EVdata(n,4)+1:EVdata(n,5),n)/0.95,
        S_EV(EVdata(n,4),n)==EVdata(n,3)+0.95*pch(EVdata(n,4),n)-pdis(EVdata(n,4),n)/0.95,
        S_EV(EVdata(n,5),n)==0.95*EVdata(n,2)];%电动汽车电量约束条件
end
ops=sdpsettings('kkt.dualbound',0);%不进行对偶边界估计
[KKTsystem,details]=kkt(C_EV,obj_EV,price_EV,ops);%将电动汽车问题转化为KKT系统
Pch=sdpvar(24,10);%储能系统充电
Pdis=sdpvar(24,10);%储能系统放电
S_ESS=sdpvar(24,10);%储能系统电量状态
Pb_DA=sdpvar(24,1);%日前购电合同
price_DA=a_data.*Pb_DA+b_data;%日前DLMP
Pb_RT=sdpvar(24,10);%实时购电量
Ps_RT=sdpvar(24,10);%实时售电量
C_price=[mean(price_EV)<=mean(price_DA),0.8*price_DA<=price_EV<=1.2*price_DA];%零售电价约束条件
C_ESS=[0<=Pch<=250,0<=Pdis<=250,200<=S_ESS<=950,
    S_ESS(1,:)==500+0.95*Pch(1,:)-Pdis(1,:)/0.95,
    S_ESS(2:24,:)==S_ESS(1:23,:)+0.95*Pch(2:24,:)-Pdis(2:24,:)/0.95,
    S_ESS(24,:)==500];%储能系统约束条件
C_CS=[0<=Pb_DA<=500,0<=Pb_RT<=500,0<=Ps_RT<=500,Pb_DA*ones(1,10)+Pb_RT+Pdis+N*pdis*ratio*ones(1,10)==Ps_RT+Pch+N*pch*ratio*ones(1,10)];%零售商约束条件
%% 求解
ops=sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.IntFeasTol',1e-9,'gurobi.MIPGap',1e-9,'gurobi.OptimalityTol',1e-9);%求解器参数,MILP问题
Constraints=[C_price,C_ESS,C_CS,KKTsystem];%总的约束条件
obj=-price_DA'*Pb_DA+sum(PDF.*sum(-(price_RT'+0.001).*Pb_RT+(price_RT'-0.001).*Ps_RT))-(details.b'*details.dual+details.f'*details.dualeq);%总的目标函数（零售商的收益）
result=optimize(Constraints,-obj,ops)%求解最大化问题
price_EV=double(price_EV);pch=double(pch);pdis=double(pdis);S_EV=double(S_EV);Pch=double(Pch);Pdis=double(Pdis);S_ESS=double(S_ESS);Pb_DA=double(Pb_DA);Pb_RT=double(Pb_RT);Ps_RT=double(Ps_RT);

