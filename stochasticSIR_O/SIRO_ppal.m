% reaction asociated (4 reactions, 5 species)
% S --> Susceptible 
% I --> Infectious (total detected WWTP)
% R --> Recovered (total)
% O --> Infected detected by Xunta
% OR -> Infected detected recovered
% N total population 
%%%%%%%% 
% S + I  --> 2I  : rc = beta*S*I/N
% I --> R        : rc = alphaI*I        %%WWTP samples
% I --> I + O    : rc = gamma*I   
% O --> OR       : rc = alphaO*I       %%Sergas samples

%% Locality selection one from [Ares,Baiona,Gondomar,Melide,Nigran]
localities = {'Ares','Baiona','Gondomar','Melide','Nigran'};
Li = 1; % 1 -> Ares; 2 -> Baiona; 3 -> Gondomar; 4 -> Melide; 5 -> Nigran 
locality = localities{Li}; % to select one of the previous localities 
 
%% Select a date with format year month day f= [yyyy mm dd] from the ones providen below
f= [2021 03 14];
%%%%%%%%%%%%%%%%%%% Posible Weeks to select %%%%%%%%%%%%%%%%%%%%%%%
%     ARES            BAIONA          GONDOMAR        MELIDE          NIGRAN
% 1	  [2020 08 29]    [2020 08 10]    [2021 01 11]    [2020 08 23]    [2020 08 15]
% 2	  [2020 09 10]    [2020 10 11]    [2021 01 18]    [2020 08 30]    [2020 09 19]
% 3	  [2020 10 17]    [2020 10 18]    [2021 01 25]    [2020 09 07]    [2020 10 12]
% 4	  [2020 10 24]    [2020 10 25]    [2021 02 01]    [2020 09 12]    [2020 10 18]
% 5	  [2020 10 31]    [2020 11 01]                    [2020 10 17]    [2020 10 25]
% 6	  [2021 01 10]    [2020 11 08]                    [2020 10 24]    [2020 11 01]
% 7	  [2021 01 16]    [2020 11 17]                    [2020 10 31]    [2020 11 21]                  
% 8	  [2021 01 24]    [2020 11 25]                    [2020 11 07]    [2020 11 28]
% 9	  [2021 01 31]    [2020 12 01]                    [2020 12 21]    [2020 12 05]
% 10  [2021 02 07]    [2020 12 05]                    [2020 12 27]    [2020 12 14]
% 11  [2021 02 28]    [2020 12 14]                    [2021 01 03]    [2020 12 21]
% 12  [2021 03 07]    [2020 12 21]                    [2021 01 10]    [2021 01 09]
% 13  [2021 03 14]    [2020 12 29]                    [2021 01 17]    [2021 01 18]
% 14  [2021 03 21]    [2021 01 04]                    [2021 01 25]    [2021 01 25]
% 15                  [2021 01 16]                    [2021 01 31]    [2021 02 01]
% 16                  [2021 01 25]                                    [2021 02 08]
% 17                  [2021 01 30]                                    [2021 02 15]
% 18                  [2021 02 06]                                    [2021 02 22]
% 19                  [2021 02 20]                                    [2021 03 01]
% 20                  [2021 03 01]                                    [2021 03 22]
% 21                  [2021 03 07]
% 22                  [2021 03 13]
% 23                  [2021 03 20]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Parameters (Pajaro et al 2022)
% Choose the ones given for each week (f) in Table 5
beta   = 0.25;  
gamma  = 0.015;
alphaI = 1/14; 
alphaO = 1/14;

%% Number of inhabitants
% Population 
PopEqui = [14493,30000,35000,14756,43000]; % population for [Ares,Baiona,Gondomar,Melide,Nigran]
hab = PopEqui(Li); 

%% Number of days to simulate
nDs = 10; % fixed in 10
Odelay = nDs; % fixed 

%% load DATA, where the number of infected persons provided by SERGAS and obtained from WWTP were saved. 
load('DATA.mat');
[~, k]=min(abs(datetime(Fecha)-datetime(f))); 
NDS = nDs -1; % since the first day data is for intinitial conditions.
B=Fecha(k:k+NDS)';
C=Fecha(k+Odelay:k+Odelay+NDS)';
if Li == 1
    %DATA for XUNTA data O observed
    RealData=I_Ares(k+Odelay:k+Odelay+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Ares(k+Odelay:k+Odelay+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
    %DATA for WWTP data I infected
    WTP_Ares(WTP_Ares == -999) = NaN;
    WTP=WTP_Ares(k:k+NDS)';
elseif Li == 2
    %DATA for XUNTA data O observed
    RealData=I_Baiona(k+Odelay:k+Odelay+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Baiona(k+Odelay:k+Odelay+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
    %DATA for WWTP data I infected
    WTP_Baiona(WTP_Baiona == -999) = NaN;
    WTP=WTP_Baiona(k:k+NDS)';
elseif Li == 3
    %DATA for XUNTA data O observed
    RealData=I_Gondomar(k+Odelay:k+Odelay+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Gondomar(k+Odelay:k+Odelay+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
    %DATA for WWTP data I infected
    WTP_Gondomar(WTP_Gondomar == -999) = NaN;
    WTP=WTP_Gondomar(k:k+NDS)';
elseif Li == 4
    %DATA for XUNTA data O observed
    RealData=I_Melide(k+Odelay:k+Odelay+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Melide(k+Odelay:k+Odelay+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
    %DATA for WWTP data I infected
    WTP_Melide(WTP_Melide == -999) = NaN;
    WTP=WTP_Melide(k:k+NDS)';
elseif Li == 5
    %DATA for XUNTA data O observed
    RealData=I_Nigran(k+Odelay:k+Odelay+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Nigran(k+Odelay:k+Odelay+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
    %DATA for WWTP data I infected
    WTP_Nigran(WTP_Nigran == -999) = NaN;
    WTP=WTP_Nigran(k:k+NDS)';
end
% check that WTP(1)>=cumI0 & WTP integers
WTP0 = WTP;
if WTP(1)<cumI0
    WTP = ceil(cumI0*WTP0/WTP0(1));
else
    WTP = ceil(WTP0);
end

%% Simulation
% Number of realizations
nsimula = 1000;
% Number of days
Tgrid=0:1:NDS; % days 
% Initial condition % [S,I,R,O,OR] 
x0 = [hab-WTP(1),WTP(1),0, cumI0,0]; 
% Maximum Infected number threshold
Imax = sum(x0(1:3));
% call the function to obtain the realizations
simulation = SIRO_ssa(alphaI,alphaO,beta,gamma,x0,Imax,Tgrid,nsimula);
II = zeros(nsimula,10);
OO = II;
newO = zeros(nsimula,10);
newO(:,1) = RealData(1);
ApproxError = zeros(nsimula,1);
ApproxErrorI = zeros(nsimula,1);

%% Results
for i=1:nsimula
    II(i,:) = simulation{i}(2,:);
    OO(i,:) = simulation{i}(4,:);
    newOT = simulation{i}(4,2:end)-simulation{i}(4,1:end-1)+...
              simulation{i}(5,2:end)-simulation{i}(5,1:end-1);
    newO(i,2:end) = newOT;
    ApproxError(i) = sum(abs(newO(i,:)-RealData));
    indO = find(WTP>=0);
    if indO > 0
        ApproxErrorI(i) = sum(abs(II(i,indO)-WTP(indO)));
    else
        ApproxErrorI(i)=0;
    end
end

acumulados =ICum0;  

figure1 = figure;
axes1 = axes('Parent',figure1,'TickLabelInterpreter','latex','FontSize',14,'Position',[0.13 0.58 0.64 0.34]);
box(axes1,'on');
hold(axes1,'on'); 
hold on
plot(B,WTP,'bo','MarkerSize',8,'LineWidth',3)
plot(B,mean(II,1),'b-','LineWidth',2)
plot(B,mean(II,1)+std(II),'b--','LineWidth',1.5)
plot(B,mean(II,1)-std(II),'b--','LineWidth',1.5)
box on
% xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
title([locality,' $I$'],'Interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
legend('Data WTP','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.81 0.77 0.17 0.15]);
% O by SERGAS 
axes2 = axes('Parent',figure1,'TickLabelInterpreter','latex','FontSize',14,'Position',[0.13 0.11 0.64 0.34]);
box(axes2,'on');
hold(axes2,'on');
plot(C,acumulados,'ks','MarkerSize',8,'LineWidth',3)
plot(C,mean(OO,1),'k-','LineWidth',2)
plot(C,mean(OO,1)+std(OO),'k--','LineWidth',1.5)
plot(C,mean(OO,1)-std(OO),'k--','LineWidth',1.5)
%xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
title([locality,' $O$'],'Interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
legend('Data O','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
legend2 = legend(axes2,'show');
set(legend2,'Position',[0.81 0.31 0.17 0.15]);
hold off
