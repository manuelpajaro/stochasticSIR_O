% We must define the number of reactions (j), the number of species (n) the
% values of the j rate constants (c0) , the initial state (n0), and the 
% state change vector (nu)

% Tgrid -> times at which we return the number of species obtained by SSA moreover
% Tgrig(1) is the initial time considered to start SSA and Tgrig(end) is
% the last time or stop condition to SSA.

function simulation = SIRv_ssa(alpha,beta0,beta1,Tchange,x0,Imax,Tgrid,nsimula)
% reaction asociated (2 reactions, 3 species)
% S --> Susceptible (S=N-I-R)
% I --> Infectious (species 1)
% R --> Recovered (species 2)
% N total population susceptible
N = sum(x0);
%%%%%%%% 
% S + I  --> 2I  : rc = beta*S*I/N
% I --> R          : rc = alpha*I
reaction_number = 2;
species_number = 3;     

nu0=zeros(reaction_number,species_number);
nu=nu0+[-1  1  0 ;
         0 -1  1]; % state change vectors dim jxn
propensity=@(x,t) (0<=t && t<Tchange)*[beta0/N*x(1)*x(2) alpha*x(2)] +...
    (t>=Tchange)*[beta1/N*x(1)*x(2) alpha*x(2)]; 

simulation=SSA_mpd_Tchange_2StopConditionALL(propensity,nu,x0,Tgrid,nsimula,Imax);
end
