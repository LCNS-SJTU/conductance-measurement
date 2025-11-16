clear;
clc;
close all;

data=load('single_dynamical_cond_estimate_tau.dat');
%data=data(:,1000:end);

for i = 0:19
    t(i+1)=i+1;
    V_w_ng(i+1,:)=data(4*i+1,:);
    V_w_g(i+1,:)=data(4*i+3,:);
    I_w_ng(i+1,:)=data(4*i+2,:);
    I_w_g(i+1,:)=data(4*i+4,:);
    diff_I(i+1,:)=I_w_ng(i+1,:)-I_w_g(i+1,:);
    Q(i+1)=sum(diff_I(i+1,1000:end));
end

ts=linspace(-100,100,length(V_w_ng(1,:)));

figure(341);
subplot(2,3,1);
plot(ts,V_w_ng');
ylabel('V without g');
xlim([-10,30]);
subplot(2,3,4);
plot(ts,V_w_g');
ylabel('V with g');
xlim([-10,30]);
ylim([-Inf,-70]);
xlabel('Time (ms)');
subplot(2,3,2);
plot(ts,I_w_ng');
ylabel('I without g');
xlim([-10,30]);
subplot(2,3,5);
plot(ts,I_w_g');
xlim([-10,30]);
xlabel('Time (ms)');
ylabel('I with g');
% Q=log(Q);
% figure(3413);
% subplot(2,1,1);
% plot(ts,diff_I');
% xlim([-10,30]);
% xlabel('Time (ms)');


subplot(2,3,[3 6]);
plot(t,Q,'o');
xlabel('Switch time (ms)');
ylabel('Mean synaptic current');

F = @(x,xdata)(x(1)*xdata+x(2)).*exp(-xdata/x(3))+x(4);
%F = @(x,xdata)(x(1).*exp(-xdata/x(2))+x(3).*exp(-xdata/x(4)));
x0 = [1 1 1 1];
x_q=Q;
x_t=t;
para1 = lsqcurvefit(F,x0,x_t,x_q);
%para2 = Diff_Evolution(x_t,x_q);

hold on;
plot(x_t,(para1(1)*x_t+para1(2)).*exp(-x_t/para1(3))+para1(4),'-r');
% plot(x_t,para2(1).*exp(-x_t/para2(2))+para2(3).*exp(-x_t/para2(4)),'-r');
