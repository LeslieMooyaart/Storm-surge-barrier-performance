clear all
close all

% 1. STORM MODELLING
% 1.1 In the analytical procedure the sea level maxima (h_sea) are the only variable of
% the storm analyzed. 

% 1.2 Exponential distribution of sea level maxima (h_sea)
h_A=2.1;%yearly maximum
h_B=0.75;%decimal height

dh_sea=h_B/750;
h_sea=0:dh_sea:h_A+10*h_B;

f_sea=log(10)/h_B*10.^(-(h_sea-h_A)/h_B);%frequency density
f_sea(h_sea<h_A)=0;

F_sea=10.^(-(h_sea-h_A)/h_B);%frequency distribution
F_sea(h_sea<h_A)=0;
E_sea=1-F_sea;%exceedance frequency

% 2. BARRIER STATE MODELLING
% 2.1 Barrier states:
%       1. Open
%       2. Failed closure
%       3. Structural failure
%       4. Succesful closure

% 2.2 Barrier state probability functions
h_c=3.0;%closure decision level [metre r.t. mean sea level]
z_str=6.6;%structural design level [metre r.t. mean sea level]
P_FC=1/100;%probability of failed closure [on demand]

for i=1:length(h_sea)
    
    if h_sea(i)>h_c
        P_1(i)=0; %open
        P_2(i)=P_FC; %failed closure
    else
        P_1(i)=1; %open
        P_2(i)=0; %failed closure
    end
    
    if h_sea(i)<z_str
        P_3(i)=0; %structural failure
    else
        P_3(i)=(1-P_2(i));%structural failure
    end
end

P_4=1-P_1-P_2-P_3; %successful closure

% 3. INNER BASIN MODELLING
% Inner water level maxima = sea level maxima

% 4. EXCEEDANCE FREQUENCY MODELLING
% Partial distributions for individual barrier modes
g_1=P_1.*f_sea; %open
g_2=P_2.*f_sea; %failed closure
g_3=P_3.*f_sea; %structural failure

F_4=sum(P_4.*f_sea)*dh_sea; %successful closure frequency
g_4=zeros(1,length(h_sea)); %length of vector similar as other partial distributions
g_4(h_c/dh_sea)=F_4/dh_sea; %successful closure

f_in=g_1+g_2+g_3+g_4; %density of hydraulic loads with barrier

for i=1:length(h_sea)
    G_1(i)=sum(g_1(i:length(h_sea))*dh_sea); %open
    G_2(i)=sum(g_2(i:length(h_sea))*dh_sea); %failed closure
    G_3(i)=sum(g_3(i:length(h_sea))*dh_sea); %structural failure    
    G_4(i)=sum(g_4(i:length(h_sea))*dh_sea); %successful closure
    F_in(i)=sum(f_in(i:length(h_sea))*dh_sea); %hydraulic loads with a barrier
end

%OUTPUT
fig1=figure(1);%exceedance frequency curve
fill_x=[F_sea+1e-10 flip(F_in)];
fill_y=[h_sea flip(h_sea)];

x1c=fill(fill_x,fill_y,[0.8 0.8 0.8],'LineStyle','none'); %surface storm surge barrier performance
hold on

G_1(G_1==0)=10^-10; %correction for asymptote for open barrier mode

x1=plot(G_1,h_sea,G_4,h_sea,G_2,h_sea,G_3,h_sea,'linewidth',2); %partial distributions barrier modes
set(x1(1), 'Color', [0 0.4470 0.7410])
set(x1(2), 'Color', [0.4660 0.6740 0.1880])
set(x1(3), 'Color', [0.8500 0.3250 0.0980])
set(x1(4), 'Color', [0.4940 0.1840 0.5560])

x1b=plot(F_in,h_sea,'k',F_sea,h_sea,'--k','linewidth',3); %hydraulic loads with and without a barrier
set(gca, 'XScale', 'log')
set ( gca, 'xdir', 'reverse')
xlabel('Exceedance frequency [per year]','fontweight','b');
ylabel('Water level maxima at Rotterdam  [metre r.t. MSL]','fontweight','b');
ylim([2.1 6.9]);
xlim([1e-7 1]);
legend([x1b(1) x1b(2) x1c(1) x1(1) x1(2) x1(3) x1(4)],'Hydraulic loads with barrier','Hydraulic loads without barrier','Performance','Open','Successful closure','Failed closure','Structural failure','Location','Northwest')

% Create textbox
annotation(fig1,'textbox',...
    [0.132142857142857 0.338095238095241 0.339285714285714 0.0928571428571452],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'String','Critical water level',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FitBoxToText','off');

% Create line
annotation(fig1,'line',[0.573214285714286 0.130357142857143],...
    [0.369047619047619 0.369047619047619],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle','-.');

% Create arrow
annotation(fig1,'arrow',[0.574999999999999 0.574999999999999],...
    [0.369047619047619 0.114285714285714],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle','-.');

fig2=figure(2);%density of inner water level maxima 
fig2.Position(3) = 800;

pos1 = [0.1 0.15 0.3 0.78];
subplot('Position',pos1);
x2a=plot(h_sea,f_sea,'--k');
hold on
x2=plot(h_sea,g_1,h_sea,g_4,h_sea,g_2,h_sea,g_3,'linewidth',2);
set(x2(1), 'Color', [0 0.4470 0.7410])
set(x2(2), 'Color', [0.4660 0.6740 0.1880])
set(x2(3), 'Color', [0.8500 0.3250 0.0980])
set(x2(4), 'Color', [0.4940 0.1840 0.5560])
xlim([2 7]);
ylim([0 3.5]);
legend([x2(1) x2(2) x2(3) x2(4) x2a(1)],'Open','Closed','Failed closure','Structural failure','Without barrier');
ylabel('Probability density [per year per meter]','fontweight','b');
xlabel('Water level maxima [meter r.t. MSL]','fontweight','b');

pos2 = [0.45 0.3 0.3 0.55];
subplot('Position',pos2);
x2b=plot(h_sea,f_sea,'--k');
hold on
x2c=plot(h_sea,g_2,h_sea,g_3,'linewidth',2);
set(x2c(1), 'Color', [0.8500 0.3250 0.0980])
set(x2c(2), 'Color', [0.4940 0.1840 0.5560])
xlim([3.6 7]);
ylim([0 3.5*10^-4]);

pos3 = [0.8 0.45 0.15 0.3];
subplot('Position',pos3);
%axes('Position',[.75 .75 .15 .15])
%box on
x2d=plot(h_sea,f_sea,'--k');
hold on
x2e=plot(h_sea,g_2,h_sea,g_3,'linewidth',2);
set(x2e(1), 'Color', [0.8500 0.3250 0.0980])
set(x2e(2), 'Color', [0.4940 0.1840 0.5560])
xlim([6 7]);
ylim([0 10^-5]);

% Create line
annotation(fig2,'line',[0.45125 0.192916666666667],...
    [0.847619047619048 0.150793650793652],'LineStyle',':');

% Create line
annotation(fig2,'line',[0.4 0.747916666666667],...
    [0.14920634920635 0.293650793650795],'LineStyle',':');

% Create line
annotation(fig2,'line',[0.80125 0.662083333333333],...
    [0.747619047619048 0.300000000000001],'LineStyle',':');

% Create line
annotation(fig2,'line',[0.949583333333333 0.750833333333333],...
    [0.446031746031746 0.3],'LineStyle',':');

% Create arrow
annotation(fig2,'arrow',[0.16 0.15875],...
    [0.858523809523813 0.933333333333337],'Color',[0.466 0.674 0.188],...
    'LineWidth',2);

% Create textarrow
annotation(fig2,'textarrow',[0.20375 0.1625],...
    [0.602380952380954 0.671428571428573],...
    'String',{'Density = \infty','Frequency =','0.06 per year'});

fig3=figure(3); % barrier state probability functions
fig3.Position(3) = 700;
Alpha=2;
Beta=2;

x3=plot(h_sea,betainv(P_4,Alpha,Beta),h_sea,betainv(P_3,Alpha,Beta),h_sea,betainv(P_2,Alpha,Beta),h_sea,betainv(P_1,Alpha,Beta),'linewidth',2);
set(x3(4), 'Color', [0 0.4470 0.7410])
set(x3(1), 'Color', [0.4660 0.6740 0.1880])
set(x3(3), 'Color', [0.8500 0.3250 0.0980])
set(x3(2), 'Color', [0.4940 0.1840 0.5560])
x3a=legend([x3(4) x3(3) x3(2) x3(1)],'open - $\overline{B}$','failed closure - $B\overline{C}$','structural failure - $BC\overline{D}$','closed - $BCD$','location','eastoutside');
set(x3a,'Interpreter','latex','fontsize',12)
ylabel('Probability [-]','fontweight','b');
xlabel('Water level maxima [metre r.t. MSL]','fontweight','b');

yticks(betainv([0 0.01 0.10 0.25 0.5 0.75 0.9 0.99, 1],Alpha,Beta))
yticklabels({'0','0.01','0.10','0.25','0.50','0.75','0.90','0.99','1'})
xlim([2 7]);

% Create textarrow
annotation(fig3,'textarrow',[0.442857142857143 0.548571428571429],...
    [0.328571428571432 0.111904761904764],...
    'String',{'Structural','failure',' level',' z_{str}'});

% Create textarrow
annotation(fig3,'textarrow',[0.241428571428571 0.185714285714286],...
    [0.206142857142858 0.111904761904762],...
    'String',{'Closure','decision','level h_c'});

Hydraulic_loads_with_barrier_table=table(F_in',h_sea');
writetable(Hydraulic_loads_with_barrier_table,"Result_storm_surge_barrier_performance_analytical.xlsx")
