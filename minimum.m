clear;
yplot2=[];
xplot2=[];
num=0;
load all_data
data_table2=data_table1;
% varTypes    = cell(1, 23);
% varTypes(1:23) = {'double'};
% varNames = {'T','dt','lmax','am','Da','D','a0','rho0','chi0','rhoc','ak','mu','ahill','gridsize','resln','exp_speed','max_rho','min_rho','min_a','zmax_zm','zmin_zm','amax','vmax'};
% data_table2 = table('Size',[1 23],'VariableTypes',varTypes,'VariableNames',varNames);
% for i=1:height(data_table1)
%     skip_val=0;
%     for j=1:height(data_table2)
%         if (ap(data_table2(j,:).Da,data_table1(i,:).Da) && ...
%             ap(data_table2(j,:).D,data_table1(i,:).D) && ap(data_table2(j,:).mu,data_table2(i,:).mu) && ...
%             ap(data_table2(j,:).chi0,data_table1(i,:).chi0) && ap(data_table2(j,:).rhoc,data_table2(i,:).rhoc) && ...
%             ap(data_table2(j,:).ak,data_table1(i,:).ak) && ap(data_table2(j,:).am,data_table1(i,:).am) &&...
%             ap(data_table2(j,:).lmax,data_table1(i,:).lmax) && data_table1(i,:).dt<data_table2(j,:).dt)
%         data_table2(i,:)=data_table1(i,:);
%         skip_val=1;
%         end
%     end
%     if skip_val==0
%         data_table2=[data_table2;data_table1(i,:)];
%     end
% end

color1='#dd614a';
color2='#f3a712';
color3='#324376';
color4='#36F1CD';
color5='#212D40';
color_green='#4f9d69';
data_table3=[];
for i=1:height(data_table2)
    if ((data_table2(i,:).chi0>=2*data_table2(i,:).D) && data_table2(i,:).rhoc>1000 && ...
            data_table2(i,:).Da>=10e-5 && (ap(data_table2(i,:).ak,data_table2(i,:).am)) &&...
            data_table2(i,:).dt<=10 && data_table2(i,:).resln>=2000 && (data_table2(i,:).min_rho>0))% &&...
          %  ap(data_table2(i,:).am,1e-3))% && ...
        %ap(data_table2(i,:).D,5/80*data_table2(i,:).Da))
        %(data_table2(i,:).chi0-data_table2(i,:).D)<=0.35*data_table2(i,:).Da)
        data_table3=[data_table3;data_table2(i,:)];
        num=num+1;
        rlc(num)=(data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da)/...
            (data_table2(i,:).chi0-data_table2(i,:).D)*(data_table2(i,:).am)/(data_table2(i,:).a0);
        yplot(num)=data_table2(i,:).min_rho;
        yplot3(num)=data_table2(i,:).max_rho;
        aplot(num)=data_table2(i,:).min_a;
        aplot2(num)=data_table2(i,:).a0/data_table2(i,:).am;
        xplot3(num)=data_table2(i,:).a0^2*data_table2(i,:).lmax/data_table2(i,:).mu...
            /data_table2(i,:).am*(data_table2(i,:).chi0-data_table2(i,:).D)/...
            (data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da);
        chiplot(num)=rlc(num)*data_table2(i,:).am*(data_table2(i,:).chi0-data_table2(i,:).D)/(data_table2(i,:).chi0);
        xplot4(num)=(data_table2(i,:).a0+(data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da)/...
            (data_table2(i,:).chi0-data_table2(i,:).D)*data_table2(i,:).am)*data_table2(i,:).lmax/data_table2(i,:).mu;
        xplot(num)=data_table2(i,:).a0*data_table2(i,:).lmax/data_table2(i,:).mu;
        xplot5(num)=(data_table2(i,:).chi0-data_table2(i,:).D)/...
            (data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da)*...
            data_table2(i,:).a0/data_table2(i,:).am;
        yplot5(num)=(yplot3(num)-yplot(num))/yplot(num);
        xplot6(num)=(data_table2(i,:).chi0)/...
            (data_table2(i,:).chi0+data_table2(i,:).Da)*...
            data_table2(i,:).a0/data_table2(i,:).am;
        xplot7(num)=(data_table2(i,:).exp_speed)^2/3600*(1/data_table2(i,:).chi0+data_table2(i,:).Da/(data_table2(i,:).chi0)^2)/...
            (data_table2(i,:).mu);
        %display((data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).Da)
        if xplot(num)>1.5*yplot(num)
            yplot2(end+1)=yplot(num);
            xplot2(end+1)=xplot(num);
            display(locn1(i))
        end
        zmax_zm(num)=data_table2(i,:).zmax_zm;
        zmin_zm(num)=-data_table2(i,:).zmin_zm;
        amax(num)=data_table2(i,:).amax;
        lambda(num)=sqrt(data_table2(i,:).lmax*data_table2(i,:).a0/data_table2(i,:).am/(data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da));
        zmax_zm_guess(num)=log(data_table2(i,:).a0/data_table2(i,:).am)/lambda(num);
        zmax_zm_guess2(num)=log(1/rlc(num))/lambda(num);
        amax_guess(num)=data_table2(i,:).a0;
        amax_guess1(num)=data_table2(i,:).a0*(data_table2(i,:).chi0)/(data_table2(i,:).chi0+data_table2(i,:).Da);
        
        amax_guess2(num)=data_table2(i,:).a0*(data_table2(i,:).chi0-data_table2(i,:).D)/(data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da);
        amax_guess3(num)=data_table2(i,:).a0*(data_table2(i,:).chi0/data_table2(i,:).D)^(-data_table2(i,:).D/(data_table2(i,:).chi0-data_table2(i,:).D));
        amax_guess4(num)=data_table2(i,:).a0*(data_table2(i,:).chi0/data_table2(i,:).D)^(-data_table2(i,:).Da/(data_table2(i,:).chi0-data_table2(i,:).D));
        vmax(num)=data_table2(i,:).vmax/data_table2(i,:).exp_speed*60;
        vmax_guess(num)=(data_table2(i,:).chi0-data_table2(i,:).D)/(data_table2(i,:).chi0-data_table2(i,:).D+data_table2(i,:).Da);
        vmax_guess(num)=rlc(num);
        ak_plot(num)=data_table2(i,:).ak/data_table2(i,:).am;
        a0(num)=data_table2(i,:).a0;
        am(num)=data_table2(i,:).am;
        chi0(num)=data_table2(i,:).chi0;
        Da(num)=data_table2(i,:).Da;
        D(num)=data_table2(i,:).D;
        
        apa(num)=data_table2(i,:).apmax/data_table2(i,:).amax;
        Nratio(num)=data_table2(i,:).Ndiff/data_table2(i,:).N0;
        chiD(num)=data_table2(i,:).D/data_table2(i,:).chi0;
        lambda2(num)=lambda(num)*(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).chi0;
        lambda3(num)=data_table2(i,:).exp_speed/data_table2(i,:).chi0/60;
        %display(data_table2(i,:).D/data_table2(i,:).Da)
        
        %sum1(num)=data_table2(i,:).Da/data_table2(i,:).chi0*data_table2(i,:).amax/(data_table2(i,:).a0-data_table2(i,:).amax);
        %sum2(num)=data_table2(i,:).mu*data_table2(i,:).Ndiff/(data_table2(i,:).a0-data_table2(i,:).amax)/data_table2(i,:).exp_speed;
        sum1(num)=data_table2(i,:).amax/data_table2(i,:).a0*(1+data_table2(i,:).Da/data_table2(i,:).chi0);
        sum2(num)=data_table2(i,:).mu*data_table2(i,:).Ndiff/data_table2(i,:).a0/(data_table2(i,:).exp_speed/60);
    end
end
%% vmax
figure(77)
clf;
semilogx(vmax_guess,1-vmax,'o')
%hold on
%loglog(vmax_guess(rlc<0.1),vmax_guess(rlc<0.1),'LineWidth',2)
xlabel('r/\lambda c')
ylabel('Simulated (c-v(z_{max}))/c')
ylim([-0.2,0.2])

%% Plots
[xplot,ord]=sort(xplot);
%% rhomin
figure(89)
clf;
loglog(xplot,xplot,'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(xplot,yplot(ord),'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('r*a_0/mu (in OD)')
ylabel('Simulated \rho_{min} (in OD)')
axis square
xlim([2e-2,5])
ylim([2e-2,5])

figure(69)
clf;
subplot(2,2,1)
cla;
loglog(xplot,xplot,'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(xplot,yplot(ord),'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('r*a_0/mu (in OD)')
ylabel('Simulated \rho_{min} (in OD)')

subplot(2,2,2)
cla;
loglog(chiplot,chiplot,'LineWidth',2);
hold on
loglog(chiplot,aplot,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('(\chi_0-D+D_a)/(\chi_0-D)*a_m^2/a_0 (in mM)')
ylabel('Simulated a_{dip} (in mM)')
ylim([5e-7,2e-4])
xlim([5e-7,2e-4])

subplot(2,2,3)
cla;
loglog(xplot3,xplot3,'LineWidth',2);
hold on
loglog(xplot3,yplot3,'o');
%loglog(xplot2,yplot2,'o');
xlabel('r*a_0^2/mu/a_m*(\chi_0-D)/(\chi_0-D+Da) (in OD)')
ylabel('Simulated \rho_{max} (in OD)')

subplot(2,2,4)
mean_fit=mean(yplot5./xplot5);
cla;
loglog(xplot5,xplot5,'LineWidth',2);
hold on
loglog(xplot5,yplot5,'o');
loglog(xplot5,xplot5/3,'LineWidth',2);
%loglog(xplot2,yplot2,'o');
xlabel('\lambda c/r')
ylabel('Simulated \rho_{max}/\rho_{min}')

figure(70)
clf;
subplot(2,1,1)
loglog(rlc,ones(size(rlc)),'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(rlc,yplot(ord)./xplot./(1-rlc),'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('r/\lambda c')
ylabel('Simulated \rho_{min}/(r*a_0/mu)')

subplot(2,1,2)
cla;
loglog(rlc,ones(size(ord)),'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(rlc,aplot./chiplot,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('r/\lambda c')
ylabel('Simulated a_{dip}/(r/\lambda_c*a_m)')
ylim([0.5,2])
xlim([5e-4,5e-1])

figure(40)
clf;
subplot(2,1,1)
loglog(ak_plot,ones(size(ak_plot)),'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(ak_plot,yplot(ord)./xplot,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('a_k/a_m')
ylabel('Simulated \rho_{min}/(r*a_0/mu)')

subplot(2,1,2)
cla;
loglog(ak_plot,ones(size(ak_plot)),'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(ak_plot,aplot./chiplot,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1);
%loglog(xplot2,yplot2,'o');
xlabel('a_k/a_m')
ylabel('Simulated a_{dip}/(r/\lambda_c*a_m)')
xlim([1e-3,1e-1])

%% rhomax/rhomin
figure(721)
mean_fit=mean(yplot5./xplot5);
clf;
loglog(xplot7,xplot7,'LineWidth',2);
hold on
loglog(xplot7,yplot5,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1,'MarkerSize',10);
%loglog(xplot6,xplot6/2,'LineWidth',2);
%loglog(xplot2,yplot2,'o');
yticks([10,100,1000])
yticklabels({'10','100','1000'})
xticks([10,100,1000])
xticklabels({'10','100','1000'})
ylim([2,2e3])
xlim([2,2e3])
axis square
set(gca,'fontsize', 32);
xlabel('\lambda c/r')
ylabel('Simulated \rho_{max}/\rho_{min}')


figure(722)
mean_fit=mean(yplot5./xplot5);
cla;
loglog(xplot6,xplot6,'LineWidth',2);
hold on
loglog(xplot6,yplot5,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1,'MarkerSize',10);
loglog(xplot6,xplot6/2,'LineWidth',2);
%loglog(xplot2,yplot2,'o');
yticks([10,100,1000])
yticklabels({'10','100','1000'})
xticks([10,100,1000])
xticklabels({'10','100','1000'})
ylim([2,2e3])
xlim([2,2e3])
axis square
set(gca,'fontsize', 32);
xlabel('\lambda c/r')
ylabel('Simulated \rho_{max}/\rho_{min}')
%%
figure(72)
subplot(2,2,2)
cla
loglog(zmax_zm_guess2,zmax_zm_guess2,'LineWidth',2)
hold on
loglog(zmax_zm_guess2,zmax_zm,'o')
xlabel('ln(\lambda c/r)/\lambda')
ylabel('Simulated z_{max}-z_m')
ylim([0.2,30])
xlim([0.2,10])

% Plot of a'/a vs lambda
subplot(2,2,3)
cla
loglog(lambda2,apa,'o')
hold on
loglog(lambda2,lambda2,'LineWidth',2)
xlabel('(\chi_0-D_\rho)\lambda/\chi_0 (analytically)')
ylabel('(da/dz)/(a(z_{max}))')

%Plot of N_max/N_0
subplot(2,2,4)
cla
loglog(xplot5,Nratio,'o')
%hold on
%loglog(chiD,chiD,'LineWidth',2)
xlabel('\lambda c/r')
ylim([0.2,0.4])
ylabel('N_1/N_0')

% Plot of a'/a vs lambda
% subplot(2,2,4)
% cla
% loglog(lambda3,apa,'o')
% hold on
% loglog(lambda3,lambda3,'LineWidth',2)
% xlabel('c/\chi_0 (numerically)')
% ylabel('(da/dz)/(a(z_{max}))')
% %xlim([0.4,40])

%% zmax_zm
figure(75)
clf
subplot(2,1,1)
zmax_zm_guess=log(xplot5)./lambda;
loglog(zmax_zm_guess2,zmax_zm_guess2,'LineWidth',2)
hold on
loglog(zmax_zm_guess2,zmax_zm,'o')
xlabel('ln(\lambda c/r)/\lambda')
ylabel('Simulated z_{max}-z_m')
legend('y=x')
ylim([0.2,6])
xlim([0.2,3])

subplot(2,1,2)
mean_fit=mean(yplot5./aplot2);
cla;
loglog(xplot5,xplot5,'LineWidth',2);
hold on
loglog(xplot5,xplot5/3,'LineWidth',2);
loglog(xplot5,yplot5,'o');
legend('y=x','y=x/3')
%loglog(xplot2,yplot2,'o');
xlim([2,2e3])
ylim([1,1e3])
xlabel('\lambda c/r')
ylabel('Simulated (\rho_{max}-\rho_{min})/\rho_{min}')

%% amax
figure(76)
clf;
subplot(2,1,1)
cla
amax_guess=(chi0)./(chi0+Da);
loglog(amax_guess,amax_guess,'LineWidth',2)
hold on
loglog(amax_guess,amax_guess/2,'LineWidth',2)
loglog(amax_guess,amax./a0,'o')
xlabel('\chi_0/(\chi_0+D_a)')
ylabel('Simulated a_{max}/a_0')
legend('y=x','y=x/2')
ylim([3e-2,2])
xlim([1e-1,1])

subplot(2,1,2)
cla
amax_guess2=a0.*(chi0-D)./(chi0-D+Da);
loglog(amax_guess2./a0,amax_guess2./a0,'LineWidth',2)
hold on
loglog(amax_guess2./a0,amax_guess2./a0/2,'LineWidth',2)
loglog(amax_guess2./a0,amax./a0,'o')
xlabel('(\chi_0-D_\rho)/(\chi_0-D_\rho+D_a)')
ylabel('Simulated a_{max}/a_0')
legend('y=x','y=x/2')
%ylim([1e-4,2])
%xlim([5e-3,1])


%% amax
figure(79)
clf;
subplot(2,1,1)
cla
loglog(amax_guess1,amax_guess1,'LineWidth',2)
hold on
loglog(amax_guess1,amax_guess1/2,'LineWidth',2)
loglog(amax_guess1,amax,'o')
xlabel('a0*(chi/D) \^ (-D/(chi-D))')
ylabel('Simulated a_{max}')
ylim([1e-3,2])
xlim([1e-2,7])

subplot(2,1,2)
cla
loglog(amax_guess4,amax_guess4,'LineWidth',2)
hold on
loglog(amax_guess4,amax,'o')
xlabel('a0*(chi/D) \^ (-Da/(chi-D))')
ylabel('Simulated a_{max}')
ylim([1e-3,2])
xlim([1e-7,7])

%% amax
figure(1)
subplot_tight(2,3,4,[0.06,0.06])
cla;
loglog(chiplot,chiplot,'LineWidth',2);
hold on
loglog(chiplot,aplot,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', 'r','MarkerSize',18);
%loglog(xplot2,yplot2,'o');
xlabel('(r/\lambda c)a_m (in nM)')
ylabel('Simulated a_{dip} (in nM)')
yticks([1,3,10,30,100]*1e-6)
yticklabels({'1','3','10','30','100'})
xticks([1,3,10,30,100]*1e-6)
xticklabels({'1','3','10','30','100'})
ylim([5e-7,2e-4])
xlim([5e-7,2e-4])
axis square
set(gca,'fontsize', 32);

subplot_tight(2,3,6,[0.06,0.06])
cla;
loglog(zmax_zm_guess2,zmax_zm_guess2,'LineWidth',2);
hold on
loglog(zmax_zm_guess2,zmin_zm,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', 'r','MarkerSize',18);
xlabel('ln(\lambda c/r)/\lambda')
ylabel('Simulated z_{min}-z_m')
yticks([0.2,0.5,1,2])
yticklabels({'0.2','0.5','1','2'})
xticks([0.2,0.5,1,2])
xticklabels({'0.2','0.5','1','2'})
ylim([0.2,2])
xlim([0.2,2])
axis square
set(gca,'fontsize', 32);

subplot_tight(2,3,5,[0.06,0.06])
cla;
loglog(xplot,xplot,'LineWidth',2);
hold on
%loglog(xplot,xplot4(ord),'LineWidth',2);
loglog(xplot,yplot(ord),'o','MarkerFaceColor', color1, 'MarkerEdgeColor', 'r','MarkerSize',18);
%loglog(xplot2,yplot2,'o');
xlabel('r*a_0/mu (in OD)')
%ylabel('Simulated \rho_{min} (in OD)')
yticks([0.03,0.1,0.3,1,3])
yticklabels({'0.03','0.1','0.3','1','3'})
xticks([0.03,0.1,0.3,1,3])
xticklabels({'0.03','0.1','0.3','1','3'})
xticks([0.03,0.3,3])
xticklabels({'0.03','0.3','3'})
axis square
xlim([2e-2,5])
ylim([2e-2,5])
set(gca,'fontsize', 32);

%%
figure(123)
mean_fit=mean(yplot5./xplot5);
cla;
loglog(xplot5,xplot5,'LineWidth',2);
hold on
loglog(xplot5,yplot5-1,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', 'r','MarkerSize',12);
%loglog(xplot5,xplot5/3,'LineWidth',2);
%loglog(xplot2,yplot2,'o');
xlabel('\lambda c/r')
ylabel('Simulated (\rho_{max}-\rho_{min})/\rho_{min}')

%%
figure(90)
clf;
%xplot5=a0./am.*(chi0-D+Da)./(chi0-D);
loglog(xplot5,sum1,'o')
hold on
loglog(xplot5,sum2,'o')
loglog(xplot5,sum1+sum2,'o')
legend('Diffusion','Consumption','Diffusion+Consumption')
xlabel('\lambda c/r')
%ylabel('Diffusion/Consumption')
ylim([0,1.1])
%title('Contribution of each term to integral of attractant equation')
%%
figure(791)
clf;
loglog(amax_guess1,amax_guess1,'LineWidth',2)
hold on
loglog(amax_guess1,amax,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', 'r','MarkerSize',12);
loglog(amax_guess1,amax_guess1/2,'LineWidth',2,'Color','k')
xlabel('a_0\chi_0/(\chi_0-D_a) (in mM)')
ylabel('Simulated a_{max} (in mM)')
ylim([1e-3,3])
xlim([1e-3,3])
axis square
set(gca,'fontsize', 32);
yticks([0.001,0.003,0.01,0.03,0.1,0.3,1,3])
yticklabels({'0.001','0.003','0.01','0.03','0.1','0.3','1','3'})
xticks([0.001,0.003,0.01,0.03,0.1,0.3,1,3])
xticklabels({'0.001','0.003','0.01','0.03','0.1','0.3','1','3'})
%%
figure(792)
clf;
zmax_zm_guess=log(xplot5)./lambda;
loglog(zmax_zm_guess2,zmax_zm_guess2,'LineWidth',2)
hold on
loglog(zmax_zm_guess2,zmax_zm,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', 'r','MarkerSize',12);
xlabel('ln(\lambda c/r)/\lambda (in mm)')
ylabel('Simulated z_{max}-z_m (in mm)')
%legend('y=x')
ylim([0.2,3])
xlim([0.2,3])
axis square
set(gca,'fontsize', 32);
yticks([0.001,0.003,0.01,0.03,0.1,0.3,1,3])
yticklabels({'0.001','0.003','0.01','0.03','0.1','0.3','1','3'})
xticks([0.001,0.003,0.01,0.03,0.1,0.3,1,3])
xticklabels({'0.001','0.003','0.01','0.03','0.1','0.3','1','3'})