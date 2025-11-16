clear;
clc;
close all;

data=load('NEURON_single_dynamical_cond_estimate_tau_only_clamp_once_E_input_dt=0.01_clampdt=0.1_len=200_V_clamp=-90_E_loc_78.dat');

dt=0.01;
delay_len=50;
clamp_dt=0.1;
for i = 0:delay_len-1
    t(i+1)=clamp_dt*(i+1);
    V_w_ng(i+1,:)=data(4*i+1,:);
    V_w_g(i+1,:)=data(4*i+3,:);
    I_w_ng(i+1,:)=data(4*i+2,:);
    I_w_g(i+1,:)=data(4*i+4,:);
    diff_I(i+1,:)=I_w_ng(i+1,:)-I_w_g(i+1,:);
    V_1(i+1,:)=zeros(1,length(V_w_ng(1,:)));
    V_1(i+1,round(5/dt)+1:round((5+clamp_dt*(i+1))/dt)+1)=V_w_g(i+1,round(5/dt)+1:round((5+clamp_dt*(i+1))/dt)+1)+70;
end

K_SS=(V_w_ng(1,end)+70)/I_w_ng(1,end);

for i = 0:delay_len-1
    Q(i+1)=sum(diff_I(i+1,:))+sum(V_1(i+1,:))/K_SS;
    Q(i+1)=Q(i+1)*dt;
end

ts=linspace(0,105,length(V_w_ng(1,:)));

figure(3413);
subplot(1,2,1);
plot(ts,diff_I');
%xlim([-10,30]);
ylabel('Synaptic current (nA)');
xlabel('Time (ms)');

subplot(1,2,2);
plot(ts,V_1');
ylabel('PSP (mV)');
%xlim([-10,30]);
xlabel('Time (ms)');



%=======================DE Method==========================================
%CASE=0;  %double exponentials
%CASE=1;  %double exponentials with fixed rise parameters
%CASE=2;  %single exponential

% paras=DE_exp_fitting(t,Q,CASE);
%  
% disp(paras);
% switch(CASE)
%         case 0
%            tau_r=1/paras(2)
%            tau_d=1/paras(4)
%         case 1
%             tau_d=1/paras(3)
%         case 2
%             tau_d=1/paras(2)
%     end
% 
% disp(K_SS);

%================second method to fit double exponentials==================
%https://math.stackexchange.com/questions/1428566/fit-sum-of-exponentials/3808325#3808325

CASE=0;  %double exponentials
%CASE=1;  %double exponentials with fixed rise parameters
%CASE=2;  %single exponential
[tau_d,tau_r,lambdas,P] = least_square_fit_multi_exponentials(t',Q',CASE);

% %================Different clamp voltage and length statistics for IF neuron=======================
% V_c=[-90,-50,-30,-10];
% E_r=[0.742,0.407,0.475,0.498];
% E_d=[2.031,1.958,1.988,1.996];
% I_r=[0.531,0.58,0.566,0.563];
% I_d=[2.005,2.015,2.013,2.012];
% figure;
% x = categorical({num2str(V_c(1)),num2str(V_c(2)),num2str(V_c(3)),num2str(V_c(4))});
% y = [E_r', I_r'];
% z=[E_d', I_d'];
% subplot(1,2,1);
% bar(x,y);
% xlabel('Clamped Voltage (mV)')
% ylabel('estimated rise time constant (ms)');
% legend('E','I')
% subplot(1,2,2);
% bar(x,z);
% xlabel('Clamped Voltage (mV)')
% ylabel('estimated decay time constant (ms)');
% legend('E','I')
% 
% 
% 
% time=[5,10,15,20,30];
% E_r=[0.5165,0.4984,0.4877,0.4805,0.4692];
% E_d=[1.9885,1.9964,1.9991,2.0004,2.0021];
% I_r=[0.5262,0.5311,0.5327,0.5333,0.5321];
% I_d=[2.0070,2.0046,2.0042,2.0041,2.0043];
% figure;
% x = categorical({num2str(time(1)),num2str(time(2)),num2str(time(3)),num2str(time(4)),num2str(time(5))});
% y = [E_r', I_r'];
% z=[E_d', I_d'];
% subplot(1,2,1);
% bar(x,y);
% xlabel('Data length (ms)')
% ylabel('estimated rise time constant (ms)');
% legend('E','I')
% subplot(1,2,2);
% bar(x,z);
% xlabel('Data length (mV)')
% ylabel('estimated decay time constant (ms)');
% legend('E','I')

% %================Different clamp voltage and length statistics for realistic neuron=======================
% 
% time=[5,10,20];
% E_d=[2.0884,2.1944,2.2718];
% E_r=[0.5109,0.5077,-0.5];
% I_d=[1.9968,1.9608,2.2531];
% I_r=[0.5126,0.5741,-0.5];
% figure;
% x = categorical({num2str(time(1)),num2str(time(2)),num2str(time(3))});
% y = [E_r', I_r'];
% z=[E_d', I_d'];
% subplot(1,2,1);
% bar(x,y);
% xlabel('Data length (ms)')
% ylabel('estimated rise time constant (ms)');
% legend('E','I')
% subplot(1,2,2);
% bar(x,z);
% xlabel('Data length (mV)')
% ylabel('estimated decay time constant (ms)');
% legend('E','I')