%% plot Everything
clear
load all_data
data_table2=data_table1;
data_table=data_table1;

%% Remove Duplicates
% varTypes    = cell(1, 23);
% varTypes(1:23) = {'double'};
% varNames = {'T','dt','lmax','am','Da','D','a0','rho0','chi0','rhoc','ak','mu','ahill','gridsize','resln','exp_speed','max_rho','min_rho','min_a','zmax_zm','zmin_zm','amax','vmax'};
% data_table2 = table('Size',[1 23],'VariableTypes',varTypes,'VariableNames',varNames);
% for i=1:height(data_table1)
%     skip_val=0;
%     for j=1:height(data_table2)
%         if (ap(data_table2(j,:).Da,data_table1(i,:).Da) && ...
%                 ap(data_table2(j,:).D,data_table1(i,:).D) && ap(data_table2(j,:).mu,data_table1(i,:).mu) && ...
%                 ap(data_table2(j,:).chi0,data_table1(i,:).chi0) && ap(data_table2(j,:).rhoc,data_table1(i,:).rhoc) && ...
%                 ap(data_table2(j,:).ak,data_table1(i,:).ak) && ap(data_table2(j,:).am,data_table1(i,:).am) &&...
%                 ap(data_table2(j,:).lmax,data_table1(i,:).lmax) && ap(data_table2(j,:).a0,data_table1(i,:).a0))
%             if data_table2(j,:).dt>=data_table1(i,:).dt
%                 data_table2(j,:)=data_table1(i,:);
%             end
%             skip_val=1;
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

%% Figure 2
% rho_height
figure(1)
clf;
% a0/am
x_plot1=[];
height_plot1=[];
locn={};
dt=[];
is=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).chi0,300e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,data_table2(i,:).am) && ap(data_table2(i,:).lmax,0.69/3600) &&...
            (data_table2(i,:).Da>1e-5) && data_table2(i,:).dt<10 && (data_table2(i,:).min_rho>0))
        height_plot1=[height_plot1,data_table2(i,:).max_rho/data_table2(i,:).min_rho-1];
        x_plot1=[x_plot1,data_table2(i,:).a0/data_table2(i,:).am];
        locn{end+1,1}=locn1{i};
        dt=[dt,data_table2(i,:).dt];
        is=[is,i];
    end
end

loglog(x_plot1,height_plot1,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
hold on

x_plot1=[];
height_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).chi0,6000e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,data_table2(i,:).am) && ap(data_table2(i,:).lmax,0.69/3600) &&...
            (data_table2(i,:).Da>1e-5) && data_table2(i,:).dt<100 && (data_table2(i,:).min_rho>0))
        height_plot1=[height_plot1,(data_table2(i,:).max_rho)/data_table2(i,:).min_rho-1];
        x_plot1=[x_plot1,data_table2(i,:).a0/data_table2(i,:).am];
    end
end

loglog(x_plot1,height_plot1,'^','MarkerFaceColor', color5, 'MarkerEdgeColor', color5,'MarkerSize',22)

x_plot1=[];
height_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ...
            ap(data_table2(i,:).D,500e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).chi0,3000e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,data_table2(i,:).am) && ap(data_table2(i,:).lmax,0.69/3600) &&...
            (data_table2(i,:).Da>1e-5) && data_table2(i,:).dt<100 && (data_table2(i,:).min_rho>0))
        height_plot1=[height_plot1,(data_table2(i,:).max_rho)/data_table2(i,:).min_rho-1];
        x_plot1=[x_plot1,data_table2(i,:).a0/data_table2(i,:).am];
    end
end

loglog(x_plot1,height_plot1,'s','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)


legend('D_\rho=50 \mu m^2/s, \chi_0=300 um^2/s, D_a=800um^2/s',...
    'D_\rho=1000 \mu m^2/s, \chi_0=6000 um^2/s, D_a=800um^2/s',...
    'D_\rho=500 \mu m^2/s, \chi_0=3000 um^2/s, D_a=800um^2/s')
ylabel('(\rho_{max}-\rho_{dip})/\rho_{dip}')
xlabel('a_0/a_m')
ylim([1,1000])
saveas(gca,'latest_figures/fig2b.png')

%% Figure 3 (r,mu,phi,Da,a0,am)
figure(2)
clf;
% r
subplot(2,2,1)

phi=5;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=linspace(0.03,3,100)/3600;
Da=800e-6;
c=D*phi*sqrt(r*a0/am/(D*phi+Da));
loglog(r*3600,c*3600,'-','Color',color1, 'LineWidth',2)
hold on
loglog(r*3600,2*sqrt(D*r)*3600,'--','Color',color1, 'LineWidth',2)
phi=5;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=linspace(.03,3,100)/3600;
Da=800e-6;
c=D*phi*sqrt(r*a0/am/(D*phi+Da));
loglog(r*3600,c*3600,'-','Color',color3, 'LineWidth',2)
loglog(r*3600,2*sqrt(D*r)*3600,'--','Color',color3, 'LineWidth',2)


x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).chi0,300e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).mu,0.77/3600) && data_table2(i,:).dt<100)
        if data_table2(i,:).exp_speed*60>7
            display(i)
        end
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).lmax];
    end
end
loglog(x_plot1*3600,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).chi0,6000e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).mu,0.77/3600) && data_table2(i,:).dt<100)
        if data_table2(i,:).exp_speed*60>7
            display(i)
        end
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).lmax];
    end
end

loglog(x_plot1*3600,vel_plot1*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
ylim([0.1,100])
xlim([0.05,1])
yticks([0.1,1,10,100])
yticklabels({'0.1','1','10','100'})
xticks([0.03, 0.1, 0.3, 1, 3])
ax = gca;
ax.FontSize = 16;
xlabel('Rate of Growth, r, in 1/h')
%ylabel('Expansion Speed (in mm/hr)')
hold on

%% mu
subplot(2,2,3)
cla;
mu=linspace(0.01,2,100)/3600;
phi=5;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=800e-6;
c=D*phi*sqrt(r*a0/am/(D*phi+Da));
semilogy(mu*3600,c*ones(size(mu))*3600,'-','Color',color1,'LineWidth',2)
hold on
semilogy(mu*3600,2*sqrt(D*r)*3600*ones(size(mu)),'--','Color',color1, 'LineWidth',2)

phi=5;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=800e-6;
c=D*phi*sqrt(r*a0/am/(D*phi+Da));
semilogy(mu*3600,c*ones(size(mu))*3600,'-','Color',color3,'LineWidth',2)
semilogy(mu*3600,2*sqrt(D*r)*3600*ones(size(mu)),'--','Color',color3, 'LineWidth',2)

x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).chi0,300e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).mu];
    end
end
semilogy(x_plot1*3600,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).chi0,6000e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        if data_table2(i,:).exp_speed*60>7
            display(i)
        end
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).mu];
    end
end
semilogy(x_plot1*3600,vel_plot1*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
yticks([0.1,1,10,100])
yticklabels({'0.1','1','10','100'})
ax = gca;
ax.FontSize = 16;
ylim([0.1,100])
%ylabel('Expansion Speed (in mm/hr)')
xlabel('Rate of uptake, \mu, in mM/OD/hr')
hold on

%% a0/am
subplot(2,2,2)
cla;
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).chi0,300e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,data_table2(i,:).am) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).a0/data_table2(i,:).am];
    end
end
mu=linspace(0,1,100);
phi=5;
D=50*1e-6;
am=1e-3;
a0=linspace(am,10,1000);
r=0.69/3600;
Da=800e-6;
c=D*phi*sqrt(r*a0/am/(D*phi+Da));
loglog(a0/am,c*3600,'-','Color',color1,'LineWidth',2)
hold on
loglog(a0/am,2*sqrt(D*r)*3600*ones(size(a0)),'--','Color',color1, 'LineWidth',2)

x_plot2=[];
vel_plot2=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).Da,800e-6) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).chi0,6000e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,data_table2(i,:).am) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        if data_table2(i,:).exp_speed*60<0.5
            display(i)
        end
        vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
        x_plot2=[x_plot2,data_table2(i,:).a0/data_table2(i,:).am];
    end
end
phi=5;
D=1000*1e-6;
am=1e-3;
a0=linspace(am,10,1000);
r=0.69/3600;
Da=800e-6;
c=D*phi*sqrt(r*a0/am/(D*phi+Da));
loglog(a0/am,c*3600,'-','Color',color3,'LineWidth',2)
ylim([0.1,100])
xlim([1,1000])
loglog(a0/am,2*sqrt(D*r)*3600*ones(size(a0)),'--','Color',color3, 'LineWidth',2)
loglog(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
loglog(x_plot2,vel_plot2*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
yticks([0.1,1,10,100])
yticklabels({'0.1','1','10','100'})
xticks([0.1,1,10,100,1000])
xticklabels({'0.1','1','10','100','1000'})
ax = gca;
ax.FontSize = 16;%ylabel('Expansion Speed (in mm/hr)')
xlabel('a_0/a_m')
hold on

%% Da
subplot(2,2,4)
cla;
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).chi0,300e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).Da];
    end
end
Da=linspace(0,2000,100)*1e-6;
phi=5;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
c=D*phi*sqrt(r*a0/am./(D*phi+Da));
semilogy(Da*1e6,c*3600,'-','Color',color1,'LineWidth',2)
hold on
semilogy(Da*1e6,2*sqrt(D*r)*3600*ones(size(Da)),'--','Color',color1, 'LineWidth',2)

x_plot2=[];
vel_plot2=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).chi0,6000e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
        x_plot2=[x_plot2,data_table2(i,:).Da];
    end
end
phi=5;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=linspace(0,2000,100)*1e-6;
c=D*phi*sqrt(r*a0/am./(D*phi+Da));
semilogy(Da*1e6,c*3600,'-','Color',color3,'LineWidth',2)
ylim([0.1,100])
xlim([20,2000])
semilogy(Da*1e6,2*sqrt(D*r)*3600*ones(size(Da)),'--','Color',color3, 'LineWidth',2)
semilogy(x_plot1*1e6,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
semilogy(x_plot2*1e6,vel_plot2*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
set(gca, 'XScale', 'log')
yticks([0.1,1,10,100])
yticklabels({'0.1','1','10','100'})
xticks([20,50,200,500,2000])
xticklabels({'20','50','200','500','2000'})
ax = gca;
ax.FontSize = 16;
%ylabel('Expansion Speed (in mm/hr)')
xlabel('Diffusion Constant of Attractant, D_a, in um^2/s')
saveas(gca,'fig3.png')
hold on

%% chi0
figure(4)
subplot(1,2,2)
cla;
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
        %x_plot1=[x_plot1,data_table2(i,:).chi0*1e6];%
    end
end
Da=800e-6;
phi=0:0.01:20;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
c1=D*phi.*sqrt(r*a0/am./(D*phi+Da));
plot(phi,c1*3600,'-','Color',color1,'LineWidth',2)
hold on
plot(phi,2*sqrt(D*r)*3600*ones(size(phi)),'--','Color',color1, 'LineWidth',2)
p=polyfit((phi(1:1000)),c1(1:1000),1);
plot(phi,p(1)*3600*phi+p(2)*3600,'--','LineWidth',10,'Color',color2)



x_plot2=[];
vel_plot2=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<100)
        vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
        x_plot2=[x_plot2,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
    end
end
phi=0:0.01:20;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=800e-6;
c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
plot(phi,c*3600,'-','Color',color3,'LineWidth',2)
ylim([0,50])
xlim([1,10])
plot(phi,2*sqrt(D*r)*3600*ones(size(phi)),'--','Color',color3, 'LineWidth',2)
p=polyfit(sqrt(phi(1000:end)),c(1000:end),1);
%p=polyfit(sqrt(x_plot2(end-10:end)),vel_plot2,1);
plot(phi,0.95*p(1)*3600*sqrt(phi)+p(2)*3600,'--','LineWidth',10, 'Color',color4)
plot(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize',10)
plot(x_plot2,vel_plot2*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize',10)


yticks([0:10:100])
%yticklabels({'0.1','1','10','100'})
ax = gca;
ax.FontSize = 16;
%xlabel('Chemotactic Sensitivity, \phi')

% axes('Position',[.18 .3 .2 .2])
%
% %ylabel('Expansion Speed (in mm/hr)')
% box on
% plot(phi(phi<1),c1(phi<1)*3600,'-','Color',color1,'LineWidth',2)
% hold on
% plot(phi(phi<1),c(phi<1)*3600,'-','Color',color3,'LineWidth',2)
% phi=-1:0.01:20;
% D=1000*1e-6;
% am=1e-3;
% a0=0.1;
% r=0.69/3600;
% Da=800e-6;
% c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
% plot(phi(phi<1),2*sqrt(1e-3*r)*3600*ones(size(phi(phi<1))),'--','Color',color3, 'LineWidth',2)
% plot(phi(phi<1),2*sqrt(50e-6*r)*3600*ones(size(phi(phi<1))),'--','Color',color1, 'LineWidth',2)
% plot(x_plot1(x_plot1<1),vel_plot1(x_plot1<1)*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
% plot(x_plot2(x_plot2<1),vel_plot2(x_plot2<1)*60,'o','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
% ax = gca;
% ax.FontSize = 16;

%% D_rho
figure(4)
subplot(1,2,1)
cla;
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(6*data_table2(i,:).D,data_table2(i,:).chi0) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).D];
    end
end
Da=800e-6;
phi=5;
D=[10:2000]*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
p=polyfit(sqrt(D(1000:end)),c(1000:end),1);
plot(D*1e6,p(1)*3600*sqrt(D)+p(2)*3600,'--','LineWidth',10, 'Color',color4)
hold on
plot(D*1e6,c*3600,'-','Color',color3,'LineWidth',2)

semilogy(x_plot1*1e6,vel_plot1*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3,'MarkerSize',10)

x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(2*data_table2(i,:).D,data_table2(i,:).chi0) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).D];
    end
end
phi=1;
D=[10:2000]*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=800e-6;
c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
p=polyfit((D(1:1000)),c(1:1000),1);
plot(D*1e6,p(1)*3600*D+p(2)*3600,'--','LineWidth',10,'Color',color2)
plot(D*1e6,c*3600,'-','Color',color1,'LineWidth',2)
ylim([0,50])
xlim([10,1200])
plot(x_plot1*1e6,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1,'MarkerSize',10)
plot(D(1:1500)*1e6,2*sqrt(D(1:1500)*r)*3600,'--','Color',color_green, 'LineWidth',3)
ylabel('Expansion Speed (in mm/hr)')
%xlabel('Motility Parameter, D_\rho')
ax = gca;
ax.FontSize = 16;
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%    'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')
saveas(gca,'fig4.png')

%% carrying capacity
figure(5)
clf
%% a0
subplot(2,2,1)
cla;
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).rhoc,10) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,300e-6))
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,data_table2(i,:).a0];
        if data_table2(i,:).a0==400
            display(locn1(i))
        end
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
a0=1e-3:1e-3:10;
r=0.69/3600;
rhoc=10;
gamma=0.26;
mu=0.77/3600;
c=D*phi*sqrt(r)./sqrt((D*phi+Da)*am./a0+r*a0./(mu*rhoc)*(D*phi)*gamma);
loglog(a0/am,c*3600,'-','Color',color1,'LineWidth',2)
hold on
loglog(a0/am,D*phi.*sqrt(r*a0/am./(D*phi+Da))*3600.*ones(size(a0)),'--','Color',color1, 'LineWidth',2)
x_plot2=[];
vel_plot2=[];
a0_max=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).rhoc,10) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,6000e-6))
        vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
        x_plot2=[x_plot2,data_table2(i,:).a0];
    end
end
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.36;
c=D*phi*sqrt(r)./sqrt((D*phi+Da)*am./a0+r*a0./(mu*rhoc)*(D*phi)*gamma);
loglog(a0/am,c*3600,'-','Color',color3,'LineWidth',2)
ylim([0.5,100])
xlim([1e-3,10]*1000)
loglog(a0/am,D*phi.*sqrt(r*a0/am./(D*phi+Da))*3600.*ones(size(a0)),'--','Color',color3, 'LineWidth',2)
loglog(x_plot1/am,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize',12)
loglog(x_plot2/am,vel_plot2*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize',12)
ylabel('Expansion Speed (in mm/hr)')
xlabel('a_0 (in mM)')
yticks([1,2,5,10,20,50])
yticklabels({'1','2','5','10','20','50'})
xticks([1,10,100,1000,1e4])
xticklabels({'0.001','0.01','0.1','1','10'})
ax = gca;
ax.FontSize = 16;
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')

%% a0_max with mu (rhoc, am, r constant)
subplot(2,2,2)
cla;
x_plot1=[];
vel_plot1=[];
a0_max=[];
mus=[];
a0s=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).rhoc,10) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).am,1e-3) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,300e-6))
        a0s=[a0s,data_table2(i,:).a0];
        match=false;
        for iii=1:length(x_plot1)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot1(iii))
                match=true;
                mus=[mus,data_table2(i,:).mu*3600];
            end
        end
        if match
            loc=find(data_table2(i,:).mu==x_plot1/(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)));
            if data_table2(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot1(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
            x_plot1=[x_plot1,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot1(a0_max==100)=[]; x_plot1(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.26;
mu=[0.01:0.01:10]/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color1,'LineWidth',2)
hold on
% subplot(2,2,4)
% cla;
% c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
% loglog(mu*rhoc/r/am,c*3600,'-','Color',color1,'LineWidth',2)
% hold on
% loglog(mu*rhoc/r/am,D*phi.*sqrt(r*a0/am./(D*phi+Da))*3600*ones(size(mu)),'--','Color',color1,'LineWidth',2)

subplot(2,2,2)
loglog(x_plot1,a0_max,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize', 12)
x_plot2=[];
vel_plot2=[];
a0_max=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).rhoc,10) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).am,1e-3) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,6000e-6))
        a0s=[a0s,data_table2(i,:).a0];
        match=false;
        for iii=1:length(x_plot2)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot2(iii))
                match=true;
                mus=[mus,data_table2(i,:).mu*3600];
            end
        end
        if match
            loc=find(data_table2(i,:).mu==x_plot2/(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)));
            if data_table2(i,:).exp_speed>vel_plot2(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot2(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
            x_plot2=[x_plot2,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot2(a0_max==100)=[];
x_plot2(a0_max==100)=[];
a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.26;
mu=[0.1:0.01:10]/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color3,'LineWidth',2)
loglog(x_plot2,a0_max,'o','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize',12)
ylabel('a_0^{max}/a_m')
xlabel('\mu \rho_c/(r a_m)')
yticks([0.05,0.1,0.2,0.5,1]*1e3)
yticklabels({'50','100','200','500','1000'})
ylim([0.05,1]*1e3)
xlim([5e3,3e4])
ax = gca;
ax.FontSize = 16;

% subplot(2,2,4)
% c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
% loglog(mu*rhoc/r/am,c*3600,'-','Color',color3,'LineWidth',2)
% loglog(mu*rhoc/r/am,D*phi.*sqrt(r*a0/am./(D*phi+Da))*3600*ones(size(mu)),'--','Color',color3,'LineWidth',2)
% loglog(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
% loglog(x_plot2,vel_plot2*60,'o','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
% ylabel('c^{max} (in mm/h)')
% xlabel('\mu \rho_c/(r a_m)')
% ylim([0.5,50])
% xlim([1e3,1e5])
% yticks([0.5,1,2,5,10,20,50])
% yticklabels({'0.5','1','2','5','10','20','50'})
% ax = gca;
% ax.FontSize = 16;
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')

subplot(2,2,2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
a0s=[];
rhocs=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).am,1e-3) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,300e-6))
        a0s=[a0s,data_table2(i,:).a0];
        rhocs=[rhocs,data_table2(i,:).rhoc];
        match=false;
        for iii=1:length(x_plot1)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot1(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot1(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
            x_plot1=[x_plot1,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot1(a0_max==100)=[]; x_plot1(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:100;
gamma=0.26;
mu=0.69/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color1,'LineWidth',2)
hold on


subplot(2,2,2)
loglog(x_plot1,a0_max,'s','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize', 12)
x_plot2=[];
vel_plot2=[];
a0_max=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).am,1e-3) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,6000e-6))
        a0s=[a0s,data_table2(i,:).a0];
        rhocs=[rhocs,data_table2(i,:).rhoc];
        match=false;
        for iii=1:length(x_plot2)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot2(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).exp_speed>vel_plot2(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot2(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
            x_plot2=[x_plot2,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot2(a0_max==100)=[]; x_plot2(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:100;
gamma=0.26;
mu=0.69/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color3,'LineWidth',2)
loglog(x_plot2,a0_max,'s','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize', 12)
ylabel('a_0^{max}/a_m')
xlabel('\mu \rho_c/(r a_m)')
yticks([0.05,0.1,0.2,0.5,1]*1e3)
yticklabels({'50','100','200','500','1000'})
ylim([0.05,1]*1e3)
xlim([5e3,3e4])
ax = gca;
ax.FontSize = 16;


%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')
subplot(2,2,2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
ams=[];
for i=1:height(data_table2)
    if (    ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).Da,800e-6) &&...
            ap(data_table2(i,:).rhoc,10) &&...
            ap(data_table2(i,:).lmax,0.69/3600) && data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,300e-6))
        a0s=[a0s,data_table2(i,:).a0];
        ams=[ams,data_table2(i,:).am];
        match=false;
        for iii=1:length(x_plot1)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot1(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot1(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
            x_plot1=[x_plot1,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot1(a0_max==100)=[]; x_plot1(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:100;
gamma=0.26;
mu=0.69/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color1,'LineWidth',2)
hold on

subplot(2,2,2)
loglog(x_plot1,a0_max,'d','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize', 12)
x_plot2=[];
vel_plot2=[];
a0_max=[];
for i=1:height(data_table2)
    if (    ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).Da,800e-6) &&...
            ap(data_table2(i,:).rhoc,10)   && ap(data_table2(i,:).lmax,0.69/3600) &&...
            data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,6000e-6))
        a0s=[a0s,data_table2(i,:).a0];
        ams=[ams,data_table2(i,:).am];
        match=false;
        for iii=1:length(x_plot2)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot2(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).exp_speed>vel_plot2(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot2(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
            x_plot2=[x_plot2,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot2(a0_max==100)=[]; x_plot2(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:100;
gamma=0.26;
mu=0.69/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color3,'LineWidth',2)
loglog(x_plot2,a0_max,'d','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize', 12)
ylabel('a_0^{max}/a_m')
xlabel('\mu \rho_c/(r a_m)')
yticks([0.05,0.1,0.2,0.5,1]*1e3)
yticklabels({'50','100','200','500','1000'})
ylim([0.1,1]*1e3)
xlim([5e3,3e4])
ax = gca;
ax.FontSize = 16;


%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')
subplot(2,2,2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
lmaxs=[];
for i=1:height(data_table2)
    if (    ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).rhoc,10) && ap(data_table2(i,:).am,1e-3) &&...
            data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,300e-6))
        a0s=[a0s,data_table2(i,:).a0];
        lmaxs=[lmaxs,data_table2(i,:).lmax];
        match=false;
        for iii=1:length(x_plot1)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot1(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot1(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
            x_plot1=[x_plot1,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot1(a0_max==100)=[]; x_plot1(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:100;
gamma=0.26;
mu=0.69/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color1,'LineWidth',2)
hold on

subplot(2,2,2)
loglog(x_plot1,a0_max,'^','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize', 12)
x_plot2=[];
vel_plot2=[];
a0_max=[];

for i=1:height(data_table2)
    if (    ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).mu,0.77/3600) && ...
            ap(data_table2(i,:).Da,800e-6) && ap(data_table2(i,:).ak,1e-3) &&...
            ap(data_table2(i,:).rhoc,10)   && ap(data_table2(i,:).am,1e-3) &&...
            data_table2(i,:).dt<100 && ap(data_table2(i,:).chi0,6000e-6))
        a0s=[a0s,data_table2(i,:).a0];
        lmaxs=[lmaxs,data_table2(i,:).lmax];
        match=false;
        for iii=1:length(x_plot2)
            if ap(data_table2(i,:).mu*(data_table2(i,:).rhoc/(data_table2(i,:).lmax*data_table2(i,:).am)),x_plot2(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).exp_speed>vel_plot2(loc)
                a0_max(loc)=data_table2(i,:).a0/data_table2(i,:).am;
                vel_plot2(loc)=data_table2(i,:).exp_speed;
            end
        else
            vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
            x_plot2=[x_plot2,(data_table2(i,:).mu*data_table2(i,:).rhoc)/(data_table2(i,:).lmax*data_table2(i,:).am)];
            a0_max=[a0_max,data_table2(i,:).a0/data_table2(i,:).am];
        end
    end
end
vel_plot2(a0_max==100)=[]; x_plot2(a0_max==100)=[]; a0_max(a0_max==100)=[];
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:100;
gamma=0.26;
mu=0.69/3600;
a0_guess=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
loglog(mu*rhoc/r/am,a0_guess/am,'-','Color',color3,'LineWidth',2)
loglog(x_plot2,a0_max,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize', 12)
ylabel('a_0^{max}/a_m')
xlabel('\mu \rho_c/(r a_m)')
yticks([0.05,0.1,0.2,0.5,1]*1e3)
yticklabels({'50','100','200','500','1000'})
xticks([3e3,5e3,1e4,3e4])
xticklabels({'3,000','5,000','10,000','30,000'})
%ylim([0.1,0.65]*1e3)
ylim([100,1e3])
xlim([2500,2e4])
ax = gca;
ax.FontSize = 16;
%legend();
saveas(gca,'fig5.png')

%%
figure(7)
subplot(1,5,[1,2,3])
cla;
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<5)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
        %display(data_table2(i,:).rhoc)
    end
end
Da=800e-6;
phi=0:0.01:20;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
c1=D*phi.*sqrt(r*a0/am./(D*phi+Da));
plot(phi,c1*3600,'-','Color',color1,'LineWidth',3)
hold on

x_plot2=[];
vel_plot2=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-2) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<100)
        vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
        x_plot2=[x_plot2,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
        %display(data_table2(i,:).rhoc)
    end
end
Da=800e-6;
[x_plot2,ord]=sort(x_plot2);
vel_plot2=vel_plot2(ord);
phi=x_plot2;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
ak=1e-2;
c1=D*phi.*sqrt(r*a0/am./(D*phi+Da));
p2=polyfit(c1(end-7:end-1),vel_plot2(end-7:end-1),1);
plot(phi,c1*60*0.63*60,'-','Color',color3,'LineWidth',3)
hold on
x_plot3=[];
vel_plot3=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-4) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<100)
        vel_plot3=[vel_plot3,data_table2(i,:).exp_speed];
        x_plot3=[x_plot3,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
        %display(data_table2(i,:).rhoc)
    end
end
[x_plot3,ord]=sort(x_plot3);
vel_plot3=vel_plot3(ord);
Da=800e-6;
phi=x_plot3;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
ak=1e-4;
c1=D*phi.*sqrt(r*a0/am./(D*phi+Da));
p1=polyfit(c1(end-7:end-1),vel_plot3(end-7:end-1),1);
plot(phi,c1*60*p1(1)*1,'-','Color',color2,'LineWidth',3)
hold on
plot(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerSize', 10)
plot(x_plot2,vel_plot2*60,'d','MarkerFaceColor', color3, 'MarkerEdgeColor', color3, 'MarkerSize', 10)
plot(x_plot3,vel_plot3*60,'s','MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerSize', 10)

yticks([0:2:20])
%yticklabels({'0.1','1','10','100'})
ax = gca;
ax.FontSize = 16;
xlim([-1,20])
%xlabel('Chemotactic Sensitivity, \phi')

% axes('Position',[.18 .3 .2 .2])
%
% %ylabel('Expansion Speed (in mm/hr)')
% box on
% plot(phi(phi<1),c1(phi<1)*3600,'-','Color',color1,'LineWidth',2)
% hold on
% plot(phi(phi<1),c(phi<1)*3600,'-','Color',color3,'LineWidth',2)
% phi=-1:0.01:20;
% D=1000*1e-6;
% am=1e-3;
% a0=0.1;
% r=0.69/3600;
% Da=800e-6;
% c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
% plot(phi(phi<1),2*sqrt(1e-3*r)*3600*ones(size(phi(phi<1))),'--','Color',color3, 'LineWidth',2)
% plot(phi(phi<1),2*sqrt(50e-6*r)*3600*ones(size(phi(phi<1))),'--','Color',color1, 'LineWidth',2)
% plot(x_plot1(x_plot1<1),vel_plot1(x_plot1<1)*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
% plot(x_plot2(x_plot2<1),vel_plot2(x_plot2<1)*60,'o','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
% ax = gca;
% ax.FontSize = 16;

subplot(1,5,[4,5])
cla;
x_plot1=[];
vel_plot1=[];
dt=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).chi0,300e-6) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt==10)
        match=false;
        for iii=1:length(x_plot1)
            if ap(data_table2(i,:).ak,x_plot1(iii))
                match=true;
                loc=iii;
            end
        end
        if match
            if data_table2(i,:).dt<dt(loc)
                x_plot1(loc)=data_table2(i,:).ak;
                vel_plot1(loc)=data_table2(i,:).exp_speed;
                dt(loc)=data_table2(i,:).dt;
            end
        else
            vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
            x_plot1=[x_plot1,data_table2(i,:).ak];
            dt=[dt,data_table2(i,:).dt];
        end
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
ak=[5e-5:1e-4:1e-1];
c1=D*phi.*sqrt(r*a0/am./(D*phi+Da));
semilogx(ak*1e3,c1*3600*ones(size(ak)),'-','Color',color1,'LineWidth',5)
hold on
semilogx(x_plot1*1e3,vel_plot1*60,'^','MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerSize', 10)
xlim([5e-5,1e-1]*1e3)
ylim([0,6])
xticks([0.01,0.1,1,10,100])
xticklabels([0.01,0.1,1,10,100])
ax = gca;
ax.FontSize = 16;
saveas(gca,'fig6.png')
%% Supplementary Inset
figure(14)
clf;
subplot(1,2,1)
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,50e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<100)
        vel_plot1=[vel_plot1,data_table2(i,:).exp_speed];
        x_plot1=[x_plot1,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
    end
end
Da=800e-6;
phi=-1:0.01:20;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
c1=D*phi.*sqrt(r*a0/am./(D*phi+Da));
semilogy(phi,c1*3600,'-','Color',color1,'LineWidth',2)
hold on
semilogy(phi,2*sqrt(D*r)*3600*ones(size(phi)),'--','Color',color1, 'LineWidth',2)

x_plot2=[];
vel_plot2=[];
for i=1:height(data_table2)
    if (ap(data_table2(i,:).mu,0.77/3600) && ap(data_table2(i,:).am,1e-3) && ...
            ap(data_table2(i,:).D,1000e-6) && ap(data_table2(i,:).a0,0.1) && ...
            ap(data_table2(i,:).Da,800e-6) && data_table2(i,:).rhoc>1000 && ...
            ap(data_table2(i,:).ak,1e-3) && ap(data_table2(i,:).lmax,0.69/3600) ...
            && data_table2(i,:).dt<100)
        vel_plot2=[vel_plot2,data_table2(i,:).exp_speed];
        x_plot2=[x_plot2,(data_table2(i,:).chi0-data_table2(i,:).D)/data_table2(i,:).D];
    end
end
phi=-1:0.01:20;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=800e-6;
c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
semilogy(phi,c*3600,'-','Color',color3,'LineWidth',2)
ylim([0.5,100])
xlim([-1,10])
semilogy(phi,2*sqrt(D*r)*3600*ones(size(phi)),'--','Color',color3, 'LineWidth',2)
semilogy(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
semilogy(x_plot2,vel_plot2*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)

yticks([0.5,2,10,100])
yticklabels({'0.5','2','10','100'})
xlabel('Chemotactic Sensitivity, \phi')
ax = gca;
ax.FontSize = 24;

ylabel('Expansion Speed (in mm/hr)')
subplot(1,2,2)
semilogy(phi(phi<1),c1(phi<1)*3600,'-','Color',color1,'LineWidth',2)
hold on
semilogy(phi(phi<1),c(phi<1)*3600,'-','Color',color3,'LineWidth',2)
phi=-1:0.01:20;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
Da=800e-6;
c=D*phi.*sqrt(r*a0/am./(D*phi+Da));
semilogy(phi(phi<1),2*sqrt(1e-3*r)*3600*ones(size(phi(phi<1))),'--','Color',color3, 'LineWidth',2)
semilogy(phi(phi<1),2*sqrt(50e-6*r)*3600*ones(size(phi(phi<1))),'--','Color',color1, 'LineWidth',2)
semilogy(x_plot1(x_plot1<1),vel_plot1(x_plot1<1)*60,'o','MarkerFaceColor', color1, 'MarkerEdgeColor', color1)
semilogy(x_plot2(x_plot2<1),vel_plot2(x_plot2<1)*60,'^','MarkerFaceColor', color3, 'MarkerEdgeColor', color3)
ylim([0.5,10])
xlabel('Chemotactic Sensitivity, \phi')
yticks([0.5,2,10])
yticklabels({'0.5','2','10'})
ylabel('Expansion Speed (in mm/hr)')
ax = gca;
ax.FontSize = 24;

%%
figure(15)
clf
% rhoc
subplot(3,2,1)
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ap(data_table(i,:).a0,0.1) && ...
            ap(data_table(i,:).D,50e-6) && ap(data_table(i,:).mu,0.77/3600) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && data_table(i,:).dt<20 && ap(data_table(i,:).chi0,300e-6))
        vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
        x_plot1=[x_plot1,data_table(i,:).rhoc];
    end
end
x_plot1(x_plot1>100)=100;
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
rhoc=1:0.1:1e2;
gamma=0.26;
mu=0.77/3600;
c=D*phi*sqrt(r)./sqrt((D*phi+Da)*am/a0+r*a0./(mu*rhoc)*(D*phi)*gamma);
plot(rhoc,c*3600,'-','Color',color1,'LineWidth',2)
hold on
plot(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
plot(rhoc,2*sqrt(D*r)*3600*ones(size(rhoc)),'--','Color',color1, 'LineWidth',2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ap(data_table(i,:).a0,0.1) &&...
            ap(data_table(i,:).D,1000e-6) && ap(data_table(i,:).mu,0.77/3600) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,6000e-6))
        vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
        x_plot1=[x_plot1,data_table(i,:).rhoc];
    end
end
x_plot1(x_plot1>100)=100;
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
a0=0.1;
r=0.69/3600;
rhoc=1:0.1:1e3;
gamma=0.36;
c=D*phi*sqrt(r)./sqrt((D*phi+Da)*am/a0+r*a0./(mu*rhoc)*(D*phi)*gamma);
plot(rhoc,c*3600,'-','Color',color3,'LineWidth',2)
ylim([0.5,50])
xlim([1,10])
plot(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
plot(rhoc,2*sqrt(D*r)*3600*ones(size(rhoc)),'--','Color',color3, 'LineWidth',2)
ylabel('Expansion Speed (in mm/hr)')
xlabel('Carrying Capacity, \rho_c')
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')

% a0
subplot(3,2,2)
x_plot1=[];
vel_plot1=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).mu,0.77/3600) && ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,50e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && data_table(i,:).dt<20 && ap(data_table(i,:).chi0,300e-6))
        vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
        x_plot1=[x_plot1,data_table(i,:).a0];
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
a0=sort(x_plot1);
r=0.69/3600;
rhoc=10;
gamma=0.26;
mu=0.77/3600;
c=D*phi*sqrt(r)./sqrt((D*phi+Da)*am./a0+r*a0./(mu*rhoc)*(D*phi)*gamma);
loglog(a0,c*3600,'-','Color',color1,'LineWidth',2)
hold on
loglog(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
loglog(a0,2*sqrt(D*r)*3600*ones(size(a0)),'--','Color',color1, 'LineWidth',2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).mu,0.77/3600) && ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,1000e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,6000e-6))
        vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
        x_plot1=[x_plot1,data_table(i,:).a0];
    end
end
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.36;
c=D*phi*sqrt(r)./sqrt((D*phi+Da)*am./a0+r*a0./(mu*rhoc)*(D*phi)*gamma);
loglog(a0,c*3600,'-','Color',color3,'LineWidth',2)
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(r*3600,c/0.01,'--','Color',color3,'LineWidth',2)
ylim([0.35,40])
xlim([1e-3,10])
loglog(x_plot1,vel_plot1*60,'o','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
loglog(a0,2*sqrt(D*r)*3600*ones(size(a0)),'--','Color',color3, 'LineWidth',2)
ylabel('Expansion Speed (in mm/hr)')
xlabel('a_0')
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')

% a0_max with mu (rhoc, am, r constant)
subplot(3,2,3)

x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,50e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,300e-6))
        if ismember(data_table(i,:).mu,x_plot1)
            loc=find(data_table(i,:).mu==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).mu];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.26;
mu=[0.1:0.01:1]/3600;
c=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
plot(mu*3600,c,'-','Color',color1,'LineWidth',2)
hold on
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(mu*3600,c*60,'--','Color',color1,'LineWidth',2)
plot(x_plot1*3600,a0_max,'o','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
plot(x_plot1*3600,vel_plot1,'s','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,1000e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,6000e-6))
        if ismember(data_table(i,:).mu,x_plot1)
            loc=find(data_table(i,:).mu-x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).mu];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.36;
mu=[0.1:0.01:1]/3600;
c=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
plot(mu*3600,c,'-','Color',color3,'LineWidth',2)
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(mu*3600,c*60,'--','Color',color3,'LineWidth',2)
ylim([0,0.6])
xlim([0.1,1])
plot(x_plot1*3600,a0_max,'o','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
plot(x_plot1*3600,vel_plot1,'s','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
ylabel('a_0^{max} (in mM) or c^{max} (in mm/min)')
xlabel('\mu')
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')
%%
% a0_max with rhoc (mu, am, r constant)
subplot(3,2,4)
%%
figure(21)
clf
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,50e-6) && ap(data_table(i,:).mu,0.77/3600) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,300e-6))
        if ismember(data_table(i,:).rhoc,x_plot1)
            loc=find(data_table(i,:).rhoc==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).rhoc];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:10;
gamma=0.26;
mu=0.77/3600;
c=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
plot(rhoc,c,'-','Color',color1,'LineWidth',2)
hold on
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(rhoc,c*60,'--','Color',color2,'LineWidth',2)
plot(x_plot1,a0_max,'o','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
plot(x_plot1,vel_plot1,'s','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,1000e-6) && ap(data_table(i,:).mu,0.77/3600) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,6000e-6))
        if ismember(data_table(i,:).rhoc,x_plot1)
            loc=find(data_table(i,:).rhoc==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).rhoc];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=0.69/3600;
rhoc=1:10;
gamma=0.36;
mu=0.77/3600;
c=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
plot(rhoc,c,'-','Color',color3,'LineWidth',2)
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(rhoc,c*60,'--','Color',color3,'LineWidth',2)
xlim([1,10])
ylim([0,0.6])
plot(x_plot1,a0_max,'o','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
plot(x_plot1,vel_plot1,'s','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
ylabel('a_0^{max} or c^{max}')
xlabel('\rho_c')
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')
% a0_max with am (rhoc, mu, r constant)
%%
subplot(3,2,5)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).mu,0.77/3600) && ...
            ap(data_table(i,:).D,50e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,data_table(i,:).am) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,300e-6))
        if ismember(data_table(i,:).am,x_plot1)
            loc=find(data_table(i,:).am==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).am];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=[0.1:0.1:10]*1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.26;
mu=0.77/3600;
c=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
semilogx(am,c,'-','Color',color1,'LineWidth',2)
hold on
c=phi*D*(mu*r*rhoc./(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
semilogx(am,c*60,'--','Color',color1,'LineWidth',2)
semilogx(x_plot1,a0_max,'o','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
semilogx(x_plot1,vel_plot1,'s','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).mu,0.77/3600) && ...
            ap(data_table(i,:).D,1000e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,data_table(i,:).am) &&...
            ap(data_table(i,:).lmax,0.69/3600) && ap(data_table(i,:).chi0,6000e-6))
        if ismember(data_table(i,:).am,x_plot1)
            loc=find(data_table(i,:).am==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).am];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=1000*1e-6;
am=[0.1:0.1:10]*1e-3;
r=0.69/3600;
rhoc=10;
gamma=0.36;
mu=0.77/3600;
c=sqrt(mu*rhoc*am/(r*gamma)*(1+Da/(D*phi)));
semilogx(am,c,'-','Color',color3,'LineWidth',2)
c=phi*D*(mu*r*rhoc./(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
semilogx(am,c*60,'--','Color',color3,'LineWidth',2)
ylim([0,1])
xlim([1e-4,1e-2])
semilogx(x_plot1,a0_max,'s','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
semilogx(x_plot1,vel_plot1,'s','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
ylabel('a_0^{max} or c^{max}')
xlabel('a_m')
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')
% a0_max with r (rhoc, am, mu constant)
subplot(3,2,6)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,50e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).mu,0.77/3600) && ap(data_table(i,:).chi0,6000e-6))
        if ismember(data_table(i,:).lmax,x_plot1)
            loc=find(data_table(i,:).lmax==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).lmax];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=50*1e-6;
am=1e-3;
r=[0.1:0.01:1]/3600;
rhoc=10;
gamma=0.26;
mu=0.77/3600;
c=sqrt(mu*rhoc*am./(r*gamma)*(1+Da/(D*phi)));
plot(r*3600,c,'-','Color',color1,'LineWidth',2)
hold on
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(r*3600,c*60,'--','Color',color1,'LineWidth',2)
plot(x_plot1*3600,a0_max,'o','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
plot(x_plot1*3600,vel_plot1,'s','MarkerFaceColor', color2, 'MarkerEdgeColor', color2)
x_plot1=[];
vel_plot1=[];
a0_max=[];
for i=1:height(data_table)
    if (ap(data_table(i,:).am,1e-3) && ...
            ap(data_table(i,:).D,1000e-6) && ap(data_table(i,:).rhoc,10) && ...
            ap(data_table(i,:).Da,800e-6) && ap(data_table(i,:).ak,1e-3) &&...
            ap(data_table(i,:).mu,0.77/3600) && ap(data_table(i,:).chi0,6000e-6))
        if ismember(data_table(i,:).lmax,x_plot1)
            loc=find(data_table(i,:).lmax==x_plot1);
            if data_table(i,:).exp_speed>vel_plot1(loc)
                a0_max(loc)=data_table(i,:).a0;
                vel_plot1(loc)=data_table(i,:).exp_speed;
            end
        else
            vel_plot1=[vel_plot1,data_table(i,:).exp_speed];
            x_plot1=[x_plot1,data_table(i,:).lmax];
            a0_max=[a0_max,data_table(i,:).a0];
        end
    end
end
Da=800e-6;
phi=5;
D=1000*1e-6;
am=1e-3;
r=[0.1:0.01:1]/3600;
rhoc=10;
gamma=0.36;
mu=0.77/3600;
c=sqrt(mu*rhoc*am./(r*gamma)*(1+Da/(D*phi)));
plot(r*3600,c,'-','Color',color3,'LineWidth',2)
c=phi*D*(mu*r*rhoc/(am*D*phi*(D*phi+Da)*gamma*4)).^(1/4);
plot(r*3600,c*60,'--','Color',color3,'LineWidth',2)
ylim([0,1])
xlim([0.1,1])
plot(x_plot1*3600,a0_max,'o','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
plot(x_plot1*3600,vel_plot1,'s','MarkerFaceColor', color4, 'MarkerEdgeColor', color4)
ylabel('a_0^{max} or c^{max}')
xlabel('r')
%legend('Analytical Expansion Speed for D_\rho=50 \mu m^2/s','Stable Fisher Speed for D_\rho=50 \mu m^2/s','Simulated Expansion Speed for D_\rho=50 \mu m^2/s',...
%   'Analytical Expansion Speed for D_\rho=1000 \mu m^2/s','Stable Fisher Speed for D_\rho=1000 \mu m^2/s','Simulated Expansion Speed for D_\rho=1000 \mu m^2/s','sd')

saveas(gca,'figS5.png')
