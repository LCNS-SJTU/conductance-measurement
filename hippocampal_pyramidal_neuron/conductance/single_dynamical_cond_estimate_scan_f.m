clear;
clc;
close all;

dt=0.1;
REV=70;

f_strength=load('E_conductance_scan_strength_active_k_h_large/strength.dat');

for trial=1:30
    xdata=0;
    xdata2=0;
    ydata=0;
    
    f_S=f_strength(trial);  %input strength
 
    tstop=50;
    data_origin=load(['E_conductance_scan_strength_active_k_h_large/',num2str(trial),'single_dynamical_cond_estimate.dat']);
    data=data_origin(:,200/dt+1:end);

    for i = 1:5
        vol(i,:)=data(1+2*i,1:end)+70;
    end

    for i = 1:5
        cur(i,:)=data(2+2*i,1:end);
    end

    func =@(x,xdata) x(1)*xdata+x(2);
    func2=@(x,xdata) x(1)*xdata;
    x0=zeros(1,2);
    x20=0;

    %===============estimate the conductance at each time step=================
    k=1;
    for i = 1:0.1/dt:tstop/dt-1
        
        for j = 1:5
            xdata(j)=vol(j,i);
            xdata2(j)=vol(j,i)-REV;
            ydata(j)=cur(j,i)-cur(j,1);
        end
        g(k,:) = lsqcurvefit(func,x0,xdata,ydata);
        g_tr(k) = lsqcurvefit(func2,x20,xdata2,ydata);

        k=k+1;
    end

    g=-g;

    kxs_kss=-g(:,1)./g(:,2)*REV;

    G_loc_rec=g(:,2)./REV./kxs_kss;
    G_loc_rec(isnan(G_loc_rec))=0;

    G_loc_tr=g_tr;
    %===============estimate the averaged conductance =========================
    k=1;
    for i = 1:0.1/dt:tstop/dt-1
        for j = 1:5
            xdata(j,k)=vol(j,i);
            ydata(j,k)=cur(j,i)-cur(j,1);
        end
        k=k+1;
    end

    m_v=mean(xdata,2);
    m_I= mean(ydata,2);
    m_g = lsqcurvefit(func,x0,m_v,m_I); 
    m_g=-m_g;

    m_kxs_kss=-m_g(:,1)./m_g(:,2)*REV;
    m_G_loc_rec(trial)=m_g(:,2)./REV./m_kxs_kss;

    %=============true local conductance===================================
    t3=linspace(0,tstop,1000);
    tr=0.5;
    td=5;
    trise=(tr*td)/(td-tr);
    factor=exp(trise/td*log(tr/td))-exp(trise/tr*log(tr/td));
    factor=1/factor;
    G_l=f_S*factor*(exp(-t3/td)-exp(-t3/tr));

    m_G_loc_tr(trial)=mean(G_loc_tr);
    true_m_G(trial)=mean(G_l);
    Error_tr(trial)=(true_m_G(trial)-m_G_loc_tr(trial))/true_m_G(trial);
    Error_m(trial)=(true_m_G(trial)-m_G_loc_rec(trial))/true_m_G(trial);
    PSP(trial)=max(abs(data(1,1:end)));
end

figure(4524352);
plot(true_m_G,m_G_loc_tr,'b.','Markersize',15);
hold on;
plot(true_m_G,m_G_loc_rec,'r.','Markersize',15);
xx=linspace(0,max(true_m_G),1000);
hold on;
plot(xx,xx,'-k');
xlabel('True mean conductance');
ylabel('Measured mean conductance');
legend('traditional method','our method');