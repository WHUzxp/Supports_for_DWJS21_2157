function result_SP=RO_SP_DC_uncertain_game(data_MP,price_RT)
%% 基础数据
N=50;%电动汽车数量
%电价数据(08:00-09:00为时段1)
%电动汽车聚类数据(08:00-09:00为时段1)
EVdata=[6,40,15,10,24,105;6,32,16,2,9,132;3,24,8,11,23,89;6,24,12,13,22,97;6,40,25,1,8,101;6,40,16,12,23,175;3,24,8,10,24,35;10,40,20,2,8,88;10,40,18,11,24,112;10,64,25,11,23,66];
%聚类得到的不同类型电动汽车的比例
ratio_initial=EVdata(:,6)/1000;%经验分布
%% 建模
pch=data_MP.pch;pdis=data_MP.pdis;price_EV=data_MP.price_EV;Pb_DA=data_MP.Pb_DA;price_DA=data_MP.price_DA;%鲁棒主问题数据
Pch=sdpvar(24,1);%储能系统充电
Pdis=sdpvar(24,1);%储能系统放电
S_ESS=sdpvar(24,1);%储能系统电量状态
Pb_RT=sdpvar(24,1);%实时购电量
Ps_RT=sdpvar(24,1);%实时售电量
ratio=sdpvar(10,1);%不同类型电动汽车的分布
C_ESS=[0<=Pch<=250,0<=Pdis<=250,200<=S_ESS<=950,
    S_ESS(1)==500+0.95*Pch(1)-Pdis(1)/0.95,
    S_ESS(2:24)==S_ESS(1:23)+0.95*Pch(2:24)-Pdis(2:24)/0.95,
    S_ESS(24)==500];%储能系统约束条件
C_CS=[0<=Pb_RT<=500,0<=Ps_RT<=500,Pb_DA+Pb_RT+Pdis+N*pdis*ratio==Ps_RT+Pch+N*pch*ratio];%零售商约束条件
obj_inner=sum(-(price_RT+0.001).*Pb_RT+(price_RT-0.001).*Ps_RT);%内层问题目标函数(最大化)
Constraints_inner=[C_ESS,C_CS];%内层问题约束条件
ops=sdpsettings('kkt.dualbound',0);%不进行对偶边界估计
[KKTsystem,details]=kkt(Constraints_inner,-obj_inner,ratio,ops);%内层问题的KKT条件
C_RO=[sum(ratio)==1,0<=ratio<=1,sum(abs(ratio-ratio_initial))<=log(20/(1-0.99))*10/2000,abs(ratio-ratio_initial)<=log(20/(1-0.99))/2000];%离散场景概率约束
%% 求解
Constraints_outer=[KKTsystem,C_RO];%外层问题约束条件
obj_outer=-price_DA'*Pb_DA+sum(-(price_RT+0.001).*Pb_RT+(price_RT-0.001).*Ps_RT)+N*price_EV'*(pch-pdis)*ratio;%外层问题目标函数(零售商的收益)
ops=sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.IntFeasTol',1e-9,'gurobi.MIPGap',1e-9,'gurobi.OptimalityTol',1e-9);%求解器参数,MILP问题
result=optimize(Constraints_outer,obj_outer,ops)%求解最小化问题
result_SP.ratio=double(ratio);result_SP.obj=double(obj_outer);
end

