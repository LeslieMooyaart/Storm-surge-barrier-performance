clear all
close all

%1. STORM MODELLING
% 1.1 In the monte carlo procedure the surge, astronomic tidal phase and river
% discharge are analyzed

% 1.2 Distributions
%Distribution maximum storm surge (dzeta_storm)
%Type Exponential
dzeta_storm_A=1.3;%
dzeta_storm_B=0.75/log(10);%

%Distribution storm duration (T_storm)
%Type lognormal
sigma_T_storm=log(1+(18.8*3600)^2/(54.3*3600)^2);%18.8;
mu_T_storm=log((54.3*3600)^2/sqrt((54.3*3600)^2+sigma_T_storm^2));

%distribution river discharge (Q);
% type lognormal
mu_Q=7.7;
sigma_Q=0.5;

%Tide
dzeta_tide=0.9;%tidal amplitude
T_tide=44712;%tidal period
b_Phi=T_tide;


% 2. BARRIER STATE MODELLING
% 2.1 Barrier states:
%       1. Open
%       2. Failed closure
%       3. Structural failure
%       4. Succesful closure

% 2.2 Barrier state probability variables
%Distribution closure decision level (h_c)
mu_h_c=3.1;
sigma_h_c=0.2;

%distribution failed closure
% Bernoulli distribution
P_FC_crit=1/100;

%Distribution structural fragility (z_str)
mu_str=6.6;
sigma_str=0.5;


% 3.1 Inner basin modelling (input)

%deterministic values
h_s=2;%start closure level
muA_mouth=3620;%river flow section
A_rhine=300e6;%Rhine basin
z_crest=5;%crest level
A_overtop=40e6;%Overflow basin

c=1.9;%overflow constant
W=360;% width barrier
g=10;% gravitational acceleration
z=3.6;%critical level

%4.1 EXCEEDANCE FREQUENCY MODELLING 
n=1e8; %number of draws in the Monte Carlo Simulation

%following to ensure all vectors are equally long
h_in_open_max=zeros(1,n);
h_in_fc_max=zeros(1,n);
h_in_str_max=zeros(1,n);
h_in_success_max=zeros(1,n);

for i=1:n
    T_storm=logninv(rand(1),mu_T_storm,sigma_T_storm);
    dzeta_storm=dzeta_storm_A-dzeta_storm_B*log(rand(1));
    h_c=norminv(rand(1),mu_h_c,sigma_h_c);
    z_str=norminv(rand(1),mu_str,sigma_str);
    P_FC=rand(1);
    Q_river=min(logninv(rand(1),mu_Q,sigma_Q),18000);%Limit river discharges to 18000 m3/s. 
    Phi=rand(1)*b_Phi;%pi*T/T_tide+1/2*pi;

    %3.2 Inner basin modelling (process)
    dt=1200;
    t=0:dt:T_storm;
    h_tide=dzeta_tide*sin(2*pi/T_tide*t-Phi);
    h_storm=dzeta_storm*cos(pi/T_storm*(t-T_storm/2)).^2;
    h_river=(8/9*Q_river/muA_mouth)^2*1/(2*g)*ones(1,length(t));
    h=h_tide+h_storm+h_river;
    h_max(i)=max(h);%
    
    %counting procedure to determine number of events over critical level
    %(z)
    v=0;
    v_str(i)=0;
    v_closed(i)=0;
    v_fc(i)=0;

    if h_max(i)>z
        v=1;
    end
    
    %2.2 barrier state and inner water levels for each monte carlo draw
    if h_max(i)<h_c %Open barrier mode
        h_in_open_max(i)=h_max(i);
        h_in_max(i)=h_in_open_max(i);
        v_open(i)=v;
    elseif P_FC<=P_FC_crit %Failed closure barrier mode
        h_in_fc_max(i)=h_max(i);
        h_in_max(i)=h_in_fc_max(i);
        v_fc(i)=v;
    elseif h_max(i)>z_str %Structural failure barrier mode
        h_in_str_max(i)=h_max(i);
        h_in_max(i)=h_in_str_max(i);
        v_str(i)=v;
    else %succesful barrier mode
        h_in_success=h;
        t_success=t;
        k = find(h>h_c,1); %first moment when closure decision level is exceeded

        Delta_h_overtop=h-z_crest;
        Delta_h_overtop=max(Delta_h_overtop,0);
        Q_overtop=c*W*Delta_h_overtop.^1.5;
        V_overtop=cumsum(Q_overtop);
        dh_overtop=V_overtop/A_overtop*dt;
        
        h_sluit=h(1:k);
        k2 = find(h_sluit<h_s,1,'last'); %moment before closure decision level when start closure level is exceeded. 
        h_initial=1;
        
        %find opening 
        h_closed=h;
        h_closed(k2:length(t))=h_initial+(t(k2:length(t))-k2*dt)*8/9*Q_river/A_rhine+dh_overtop(k2:length(t));
        k3 = find(h-h_closed<0,1);
               
        %inner water level
        h_in_success(k2:k3)=h_closed(k2:k3);
        h_in_success_max(i)=min(max(h_in_success),h_max(i));%Open with very extreme river discharge
        h_in_max(i)=h_in_success_max(i);
        
        
    end

    F(i)=i/n;     
end

Sorted_h_without=sort(h_max,'descend');%without barrier
Sorted_h_in_max=sort(h_in_max,'descend');
Sorted_h_in_open_max=sort(h_in_open_max,'descend');
Sorted_h_in_success_max=sort(h_in_success_max,'descend');
Sorted_h_in_fc_max=sort(h_in_fc_max,'descend');
Sorted_h_in_str_max=sort(h_in_str_max,'descend');

F_crit_open=sum(v_open)/n;
F_crit_fc=sum(v_fc)/n;
F_crit_str=sum(v_str)/n;
F_crit_closed=sum(v_closed)/n;



%% Only change output figure 1

%Probability functions plot
H_2=2:0.1:8;%variable used to explain the input of barrier state probability functions
P_CDL=1-normcdf(H_2,mu_h_c,sigma_h_c);%Open
P_FC2=P_FC_crit.*(1-P_CDL);%Failed closure
P_str=(1-P_CDL-P_FC2).*normcdf(H_2,mu_str,sigma_str);%Structural failure
P_success=1-P_CDL-P_FC2-P_str;

Alpha=2;
Beta=2;

fig1=figure(1);
fig1.Position(3) = 700;
x1=plot(H_2,betainv(P_CDL,Alpha,Beta),H_2,betainv(P_FC2,Alpha,Beta),H_2,betainv(P_str,Alpha,Beta),H_2,betainv(P_success,Alpha,Beta),'linewidth',2);
set(x1(1), 'Color', [0 0.4470 0.7410])
set(x1(4), 'Color', [0.4660 0.6740 0.1880])
set(x1(2), 'Color', [0.8500 0.3250 0.0980])
set(x1(3), 'Color', [0.4940 0.1840 0.5560])
x1a=legend('open - $\overline{B}$','failed closure - $B\overline{C}$','structural failure - $BC\overline{D}$','closed - $BCD$','location','eastoutside');
set(x1a,'Interpreter','latex','fontsize',12)

ylabel('Probability [-]','fontweight','b')
xlabel('Sea water level maximum [m r.t. MSL]','fontweight','b')
yticks(betainv([0 0.01 0.10 0.25 0.5 0.75 0.9 0.99, 1],Alpha,Beta))
yticklabels({'0','0.01','0.10','0.25','0.50','0.75','0.90','0.99','1'})


%% Only change figure 2

A=readtable("Result_storm_surge_barrier_performance_analytical.xlsx");
F_in=A.Var1;
h_sea=A.Var2;

fig2=figure(2);

F_red=[F(1:1000) F(1001:1000:n)];
Sorted_H_red=[Sorted_h_without(1:1000) Sorted_h_without(1001:1000:n)];
Sorted_h_in_max_red=[Sorted_h_in_max(1:1000) Sorted_h_in_max(1001:1000:n)];

fill_x=[F_red flip(F_red)];
fill_y=[Sorted_H_red flip(Sorted_h_in_max_red)];

x3c=fill(fill_x,fill_y,[0.8 0.8 0.8],'LineStyle','none');
hold on
x3a=plot(F_in,h_sea,':',F,Sorted_h_in_open_max,F,Sorted_h_in_success_max,F,Sorted_h_in_fc_max,F,Sorted_h_in_str_max,'linewidth',2);
x3=plot(F,Sorted_h_without,'--k',F,Sorted_h_in_max,'k','linewidth',3);
set(gca, 'XScale', 'log', 'xdir', 'reverse')
set(x3a(2), 'Color', [0 0.4470 0.7410])
set(x3a(3), 'Color', [0.4660 0.6740 0.1880])
set(x3a(4), 'Color', [0.8500 0.3250 0.0980])
set(x3a(5), 'Color', [0.4940 0.1840 0.5560])
set(x3a(1), 'Color', [0 0 0])
legend([x3(1) x3(2) x3a(1) x3c(1) x3a(2) x3a(3) x3a(4) x3a(5)],'without barrier','with barrier','with barrier (analytical)','performance','open','closed','failed closure','structural failure','location','northwest')
ylim([2 7]);
xlim([1e-7 1]);
xlabel('Exceedance frequencies [per year]','fontweight','b')
ylabel('Water level maximum at Rotterdam [m r.t. MSL]','fontweight','b')

% Create line
annotation(fig2,'line',[0.132142857142857 0.748214285714285],...
    [0.368047619047619 0.369047619047619],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle','-.');

% Create arrow
annotation(fig2,'arrow',[0.557142857142856 0.557142857142856],...
    [0.368047619047619 0.10952380952381],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle','-.');

% Create textbox
annotation(fig2,'textbox',...
    [0.127785714285714 0.376190476190477 0.256142857142857 0.0523809523809523],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'String','Critical water level',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FitBoxToText','off');

% Create arrow
annotation(fig2,'arrow',[0.748214285714286 0.748214285714285],...
    [0.371428571428571 0.111904761904762],...
    'HeadSize',8,'LineStyle','--','Color', [0.4940 0.1840 0.5560]);

% Create arrow
annotation(fig2,'arrow',[0.68392857142857 0.68392857142857],...
    [0.36904761904762 0.10952380952381],...
    'HeadSize',8,'LineStyle','--','Color',[0 0.4470 0.7410]);

% Create arrow
annotation(fig2,'arrow',[0.64642857142857 0.64642857142857],...
    [0.36904761904762 0.10952380952381],...
    'HeadSize',8,'LineStyle','--','Color', [0.4660 0.6740 0.1880]);

% Create arrow
annotation(fig2,'arrow',[0.574999999999998 0.574999999999998],...
    [0.36904761904762 0.109523809523811],...
    'HeadSize',8,'LineStyle','--','Color', [0.8500 0.3250 0.0980] );

%% Only change figure 3

A1=nonzeros(h_in_open_max);
A2=nonzeros(h_in_success_max);
A3=nonzeros(h_in_fc_max);
A4=nonzeros(h_in_str_max);

binrange1=2.1:0.1:7;%create binrange

countsa=histc(A1,binrange1);
countsb=histc(A2,binrange1);
countsc=histc(A3,binrange1);
countsd=histc(A4,binrange1);

bar_open=countsa+countsb+countsc+countsd;
bar_closed=countsb+countsc+countsd;
bar_failed=countsc+countsd;
bar_sf=countsd;


fig3=figure(3);
fig3.Position(3) = 800;

pos1 = [0.1 0.15 0.3 0.78];
subplot('Position',pos1);
x4a=bar(binrange1,bar_open,'linestyle','none');
set(x4a(1), 'FaceColor', [0 0.4470 0.7410])
hold on
x4b=bar(binrange1,bar_closed,'linestyle','none');
set(x4b(1), 'FaceColor', [0.4660 0.6740 0.1880])
x4c=bar(binrange1,bar_failed,'linestyle','none');
set(x4c(1), 'FaceColor', [0.8500 0.3250 0.0980])
x4d=bar(binrange1,bar_sf,'linestyle','none');
set(x4d(1), 'FaceColor', [0.4940 0.1840 0.5560])
xlim([1 7])
legend([x4a x4b x4c x4d],'open','closed','failed closure','structural failure','location','northeast')
ylabel('Number of samples [per year per meter]','fontweight','b');
xlabel('Water level maxima [meter r.t. MSL]','fontweight','b');
hold off

pos2 = [0.45 0.3 0.3 0.55];
subplot('Position',pos2);
x4e=bar(binrange1,bar_open,'linestyle','none');
set(x4e(1), 'FaceColor', [0 0.4470 0.7410])
hold on
x4f=bar(binrange1,bar_closed,'linestyle','none');
set(x4f(1), 'FaceColor', [0.4660 0.6740 0.1880])
x4g=bar(binrange1,bar_failed,'linestyle','none');
set(x4g(1), 'FaceColor', [0.8500 0.3250 0.0980])
x4h=bar(binrange1,bar_sf,'linestyle','none');
set(x4h(1), 'FaceColor', [0.4940 0.1840 0.5560])
xlim([z 7])
ylim([0 500*n/1e7])
%ylabel('Number of samples [per year per meter]','fontweight','b');
%xlabel('Water level maxima [meter r.t. MSL]','fontweight','b');

pos3 = [0.8 0.45 0.15 0.3];
subplot('Position',pos3);
x4i=bar(binrange1,bar_open,'linestyle','none');
set(x4i(1), 'FaceColor', [0 0.4470 0.7410])
hold on
x4j=bar(binrange1,bar_closed,'linestyle','none');
set(x4j(1), 'FaceColor', [0.4660 0.6740 0.1880])
x4k=bar(binrange1,bar_failed,'linestyle','none');
set(x4k(1), 'FaceColor', [0.8500 0.3250 0.0980])
x4l=bar(binrange1,bar_sf,'linestyle','none');
set(x4l(1), 'FaceColor', [0.4940 0.1840 0.5560])
xlim([5 7])
ylim([0 10*n/1e7])
% Create line
annotation(fig3,'line',[0.22875 0.450416666666667],...
    [0.154761904761905 0.843650793650794],'LineStyle',':');

% Create line
annotation(fig3,'line',[0.4 0.74875],...
    [0.152380952380952 0.302380952380952],'LineStyle',':');

% Create line
annotation(fig3,'line',[0.74875 0.95],...
    [0.302380952380952 0.452380952380952],'LineStyle',':');

% Create line
annotation(fig3,'line',[0.57125 0.80125],...
    [0.311904761904763 0.747619047619048],'LineStyle',':');




