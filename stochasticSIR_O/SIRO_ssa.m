% We must define the number of reactions (j), the number of species (n) the
% values of the j rate constants (c0) , the initial state (n0), and the 
% state change vector (nu)

% Tgrid -> times at which we return the number of species obtained by SSA moreover
% Tgrig(1) is the initial time considered to start SSA and Tgrig(end) is
% the last time or stop condition to SSA.

function simulation = SIRO_ssa(alphaI,alphaO,beta,gamma,x0,Imax,Tgrid,nsimula)
% reaction asociated (4 reactions, 4 species)
% S --> Susceptible (S=N-I-R)
% I --> Infectious (species 1)
% R --> Recovered (species 2)
% O --> Infected detedted by Xunta
% N total population susceptible
N = sum(x0(1:3));
%%%%%%%% 
% S + I  --> 2I  : rc = beta*S*I/N
% I --> R        : rc = alphaI*I
% I --> I + O    : rc = gamma*I  
% O --> 0        : rc = alphaO*O
reaction_number = 4;
species_number = 5;     

c0=zeros(1,reaction_number);
c0=c0+[beta/N alphaI gamma alphaO];
nu0=zeros(reaction_number,species_number);
nu=nu0+[-1  1  0  0 0;
         0 -1  1  0 0;
         0  0  0  1 0;
         0  0  0 -1 1]; % state change vectors dim jxn
propensity=@(x) [c0(1)*x(1)*x(2), c0(2)*x(2), c0(3)*x(2), c0(2)*x(4)]; 

simulation=SSA_mpd_2StopConditionALL(propensity,nu,x0,Tgrid,nsimula,Imax);
end
