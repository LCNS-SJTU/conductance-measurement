clear;
clc;
close all;

dt=0.1;    %simulation time used in NEURON

REV=70;    %reversal potential; 70mV for E and -10mV for I
f_S=1e-3;  %input strength

tstop=100;
data_origin=load('single_dynamical_cond_estimate.dat');
data=data_origin(:,199/dt+1:end);
vol=zeros(5, 1011);
cur=zeros(5, 1011);

data(1:2:end,:)=data(1:2:end,:)+70;
figure(1);
for i = 1:5
    vol(i,:)=data(1+2*i,1:end);
%     plot(vol(i,:));
%     hold on;
end
% figure(2);
for i = 1:5 
    cur(i,:)=data(2+2*i,1:end);
%     plot(cur(i,:))
%     hold on
end

% func =@(x,xdata) x(1)*xdata+x(2);
% func2=@(x,xdata) x(1)*xdata;
% x0=zeros(1,2);

% %===============estimate the conductance at each time step=================
% k=1;
% for i = 1:0.1/dt:tstop/dt-1
%     for j = 1:5
%         xdata(j)=vol(j,i);
%         xdata2(j)=vol(j,i)-REV;
%         ydata(j)=cur(j,1)-cur(j,i);
%     end
%     g(k,:) = lsqcurvefit(func,x0,xdata,ydata);     %k and b measured by our method
%     g_tr(k) = lsqcurvefit(func2,x0(1),xdata2,ydata); %k and b measured by the traditional method
%     
%     %==========plot the IV relation I=g(V-Erev) at t=15ms==================
%     if k==10/dt
%         figure(7442342);
%         xxxxx=xdata2;
%         yyyyy=ydata;
%         plot(xdata2,ydata,'o');
%         hold on;
%         xx=linspace(min(xdata2),max(xdata2),1000);
%         yy=g_tr(k)*xx;
%         plot(xx,yy,'-r');
%         xlabel('Vol-RevE (mV)');
%         ylabel('synaptic current (nA)');
%         title('t=10 ms');
%     end
%     
%     k=k+1;
% end
% 
% kes_kss=-g(:,1)./g(:,2)*REV;
% 
% G_loc_rec=g(:,2)./REV./kes_kss;  %recovered local GE at each time point 
% G_loc_rec(isnan(G_loc_rec))=0;
% 
% G_loc_tr=-g_tr; %recovered local GE at each time point using the traditional method

%===============estimate the averaged conductance =========================

xdata=zeros(5,999);
ydata=zeros(5,999);
k=1;
for i = 1:0.1/dt:tstop/dt-1
    for j = 1:5
        xdata(j,k)=vol(j,i);
        ydata(j,k)=cur(j,1)-cur(j,i);
    end
    k=k+1;
end

m_v=mean(xdata,2);
m_I= mean(ydata,2);
% m_g = lsqcurvefit(func,x0,m_v,m_I);
% 
% figure(33241);
% plot(m_v,m_I,'o');
% hold on;
% xx=linspace(min(xdata(:,1)),max(xdata(:,1)),1000);
% yy=m_g(1)*xx+m_g(2);
% plot(xx,yy,'-r');
% xlabel('holding votlage (mV)');
% ylabel('mean synaptic current (nA)');
% 
% m_kxs_kss=-m_g(:,1)./m_g(:,2)*REV;
% m_G_loc_rec=m_g(:,2)./REV./m_kxs_kss  %mean local conductance measured by our method 
% 
% %=============plot results===================================
% figure(1);
% t2=linspace(-15,tstop,length(G_loc_rec));
% plot(t2,G_loc_rec,'-b');
% hold on;
% plot(t2,G_loc_tr,'--b');
% t3=linspace(0,tstop,1000);
% 
% %==================true local condutance==========================
% tr=0.5;
% td=5;
% trise=(tr*td)/(td-tr);
% factor=exp(trise/td*log(tr/td))-exp(trise/tr*log(tr/td));
% factor=1/factor;
% G_l=f_S*factor*(exp(-t3/td)-exp(-t3/tr));
% plot(t3,G_l,'-r');
% ylim([0,Inf]);
% legend('our method','traditional method','true conductance');
% xlabel('Time (ms)');
% ylabel('Conductance (nS)');
% xlim([0,tstop]);
% %=================================================================
% 
% figure(2)
% plot(t2,kes_kss);
% ylim([0,1]);
% xlabel('Time (ms)');
% ylabel('Factor');
% 
% m_G_loc_tr=mean(G_loc_tr)        %mean local condutance measured by the traditional method
% true_m_G=mean(G_l)               %mean true local conductance 
% Error_tr=(true_m_G-m_G_loc_tr)/true_m_G %error of mean local condutance measured by the traditional method
% Error_m=(true_m_G-m_G_loc_rec)/true_m_G  %error of mean true local conductance
% 
% t4=linspace(-15,tstop,length(vol(1,:)));
% figure(3)
% subplot(1,3,1);
% plot(t4,vol,'Linewidth',2);
% xlim([-15,tstop]);
% xlabel('Time (ms)');
% ylabel('Voltage (mV)');
% title('clamp voltage');
% subplot(1,3,2);
% plot(t4,cur,'Linewidth',2)
% xlim([-15,tstop]);
% xlabel('Time (ms)');
% ylabel('Current (nA)');
% title('clamp current');
% subplot(1,3,3);
% plot(t4,data(1,1:end),'-r','Linewidth',2)
% xlim([-15,tstop]);
% xlabel('Time (ms)');
% ylabel('Voltage (mV)');
% title('PSP');