clear;
clc;
close all;

data=load('volt_spatial_distribution_vclamp.dat');
dist=load('distance.dat');
dist(1)=0;
data=data(:,end-100)+70;
for i =1:10
    vol(i)=-110+i*10;
    v_x(i,:)=data((1:18)+(i-1)*18);
    figure(21341);
    subplot(1,2,1);
    plot(dist,v_x(i,:),'-');
    hold on;
    subplot(1,2,2);
    plot(dist,v_x(i,:)/v_x(i,1),'-');
    hold on;
end
figure(21341);
subplot(1,2,1);
legend('-100','-90','-80','-70','-60','-50','-40','-30','-20','-10');
xlabel('distance (\mu m)');
ylabel('Voltage (mV)');
subplot(1,2,2);
legend('-100','-90','-80','-70','-60','-50','-40','-30','-20','-10');
xlabel('distance (\mu m)');
ylabel('K_{ES}/K_{SS}');
