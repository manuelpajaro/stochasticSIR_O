% reaction asociated (6 reactions, 3 species)
% S --> Susceptible 
% I --> Infectious 
% R --> Recovered 
% N total population
%%%%%%%% 
% S + I  --> 2I  : rc = beta*S*I/N
% 2S + I --> 3I  : rc = beta*S*(S-1)*I/N^2
% 3S + I --> 4I  : rc = beta*S*(S-1)*(S-2)*I/N^3
% 4S + I --> 5I  : rc = beta*S*(S-1)*(S-2)*(S-3)*I/N^4
% 5S + I --> 6I  : rc = beta*S*(S-1)*(S-2)*(S-3)*(S-4)*I/N^5
% I --> R          : rc = alpha*I

%% Locality selection one from [Ares,Baiona,Gondomar,Melide,Nigran]
localities = {'Ares','Baiona','Gondomar','Melide','Nigran'};
Li = 4; % 1 -> Ares; 2 -> Baiona; 3 -> Gondomar; 4 -> Melide; 5 -> Nigran 
locality = localities{Li}; % to select one of the previous localities 

%% Select a date with format year month day f= [yyyy mm dd] from the ones providen below
f= [2020 10 25];
%%%%%%%%%%%%%%%%%%% Posible Weeks to select %%%%%%%%%%%%%%%%%%%%%%%
%  ARES            BAIONA          GONDOMAR        MELIDE          NIGRAN
% [2020 08 30]    [2021 03 14]    [2020 08 07]    [2020 03 20]    [2020 03 12]
% [2020 10 18]                    [2020 10 11]    [2020 08 09]    [2020 08 02]
% [2021 01 17]                    [2020 11 22]    [2020 08 16]    [2021 03 14]
% [2021 04 11]                    [2021 03 14]    [2020 10 25]    
%                                                 [2021 01 03]

%% Model Parameters (Pajaro et al Environ. Modell. Softw. 2022)
% Choose the ones given for each week (f) in Table 4
beta  = 0.02; % beta in {0.02,0.03}
alpha = 1/14; % fixed

%% Number of inhabitants 
Population = [14493,30000,35000,14756,43000]; % population in [Ares,Baiona,Gondomar,Melide,Nigran]
hab = Population(Li); 

%% Number of days to simulate
nDs = 7; % fixed in 7

%% load DATA, where the number of infected persons provided by SERGAS were saved. 
load('DATA.mat');
NDS = nDs -1; % since the first day data is for intinitial conditions.
[~, k]=min(abs(datetime(Fecha)-datetime(f))); %being datetime(year,month,day)
A=Fecha(k:k+NDS)';
if Li == 1
    RealData=I_Ares(k:k+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Ares(k:k+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
elseif Li == 2
    RealData=I_Baiona(k:k+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Baiona(k:k+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
elseif Li == 3
    RealData=I_Gondomar(k:k+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Gondomar(k:k+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
elseif Li == 4
    RealData=I_Melide(k:k+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Melide(k:k+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
elseif Li == 5
    RealData=I_Nigran(k:k+NDS)';%Real Data Infected detecte by Xunta
    ICum0=ICum14_Nigran(k:k+NDS)'; %cummulative ifected 14 days
    cumI0 =ICum0(1); 
end

%% Simulation
% Number of realizations
nsimula = 1000;
% Number of days
Tgrid=0:1:NDS; % days
% Initial condition 
x0 = [hab-cumI0,cumI0,0];  
Imax = 500;% Maximum Infected number threshold
% call the function to obtain the realizations
restriction = 0; %DO NOT MODIFY
simulation = SIRj_ssa(alpha,beta,restriction,x0,Imax,Tgrid,nsimula);

%% plot results
sizeL = 14; % size for figures axis,labels,titles
II = zeros(nsimula,nDs);
Inew = zeros(nsimula,nDs);
Inew(:,1) = RealData(1);
ApproxError = zeros(nsimula,1);

for i=1:nsimula
    II(i,:)=simulation{i}(2,:);
    InewT = simulation{i}(1,1:end-1)-simulation{i}(1,2:end);
    Inew(i,2:end) = InewT;
    ApproxError(i) = sum(abs(Inew(i,:)-RealData));
end

figure1=figure;
hold on
%axes1 = axes('Parent',figure1,'TickLabelInterpreter','latex','FontSize',14,'Position',[0.13 0.13 0.7 0.8]);
%hold(axes1,'on');
%plot(A,ICum0,'k-','LineWidth',2)
plot(A,ICum0,'ks','MarkerSize',8,'LineWidth',3)
plot(A,mean(II,1),'b-','LineWidth',2)
plot(A,mean(II,1)+std(II),'b--','LineWidth',1.5)
plot(A,mean(II,1)-std(II),'b--','LineWidth',1.5)
%xlabel('Day','Fontsize',sizeL,'Interpreter','latex')
ylabel(' Infected','Fontsize',sizeL,'Interpreter','latex')
title([locality,', ',num2str(f(1))],'Interpreter','latex')
legend('Real','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
%legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.84 0.4 0.15 0.50]);
set(gca,'FontSize',sizeL,'TickLabelInterpreter','latex')
hold off

%%number of realizations that are similar to the given data
ind0 = find(ApproxError==0);
ind1 = find(ApproxError==1);
ind2 = find(ApproxError==2);
ind3 = find(ApproxError==3);

%%plot the best realization
bestind = find(ApproxError==min(ApproxError));
figure
hold on
plot(A,RealData,'ks','MarkerSize',8,'LineWidth',3)
plot(A,Inew(bestind(1),:),'rs','MarkerSize',3,'LineWidth',3)
%xlabel('Day','Fontsize',sizeL,'Interpreter','latex')
ylabel('New Infected','Fontsize',sizeL,'Interpreter','latex')
legend('Real','SSA')
title(['SIRj, ',locality,', ',num2str(f(1))],'Interpreter','latex')
set(gca,'FontSize',sizeL,'TickLabelInterpreter','latex')
ylim([0 max([RealData,Inew(bestind(1),:)])+1])
hold off

%% histograms
hist2plot = 1; % 0 -> new infected cases; 1 -> cumulate infected cases 
if hist2plot == 0
    histvar = Inew;
    histdat = RealData;
    xnamevar = '$I_{New}$';
    ynamevar = '$P(I_{New})$';
else
    histvar = II;
    histdat = ICum0;
    xnamevar = '$I_{Cum}$';
    ynamevar = '$P(I_{Cum})$';
end
Imin=min(min(histvar));
Imax=max(max(histvar));
figure('Renderer', 'painters', 'Position', [100 100 800 500])
for i=2:nDs
    subplot(2,3,i-1)
    hold on
    hhI=histogram(histvar(:,i),'Normalization','pdf');
    plot([histdat(i),histdat(i)],[0,1],'k:','LineWidth',4)
    ylim([0,1.3*max(hhI.Values)])
    xlim([Imin,Imax])
    xlabel(xnamevar,'Fontsize',sizeL,'Interpreter','latex')
    ylabel(ynamevar,'Fontsize',sizeL,'Interpreter','latex')
    title([locality,', ',datestr(A(i))],'Interpreter','latex')
    set(gca,'FontSize',sizeL,'TickLabelInterpreter','latex')
    hold off
end
