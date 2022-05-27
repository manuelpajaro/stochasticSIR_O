% We must define the number of reactions (j), the number of species (n) the
% values of the j rate constants (c0) , the initial state (n0), and the 
% state change vector (nu)

% Tgrid -> times at which we return the number of species obtained by SSA moreover
% Tgrig(1) is the initial time considered to start SSA and Tgrig(end) is
% the last time or stop condition to SSA.

function simulation = SIRj_ssa(alpha,beta,restriction,x0,Imax,Tgrid,nsimula)
% reaction asociated (2 reactions, 3 species)
% S --> Susceptible (S=N-I-R)
% I --> Infectious (species 1)
% R --> Recovered (species 2)
% N total population susceptible
N = sum(x0);
%%%%%%%% 
% S + I  --> 2I  : rc = beta*S*I/N
% 2S + I --> 3I  : rc = beta*S*(S-1)*I/N^2
% 3S + I --> 4I  : rc = beta*S*(S-1)*(S-2)*I/N^3
% 4S + I --> 5I  : rc = beta*S*(S-1)*(S-2)*(S-3)*I/N^4
% 5S + I --> 6I  : rc = beta*S*(S-1)*(S-2)*(S-3)*(S-4)*I/N^5
% I --> R          : rc = alpha*I
reaction_number = 6;
species_number = 3;     

% exponential rate
if restriction == 1
    R0 = 2; % 1.5<R0<3.5
    yy=1:5;
    fyy = 1/R0*exp(-yy/R0);
else
    fyy = ones(1,5);
end

c0=zeros(1,reaction_number);
c0=c0+[beta/N*fyy(1) beta/N^2*fyy(2) beta/N^3*fyy(3) beta/N^4*fyy(4) beta/N^5*fyy(5) alpha];
nu0=zeros(reaction_number,species_number);
nu=nu0+[-1  1  0 ;
        -2  2  0 ;
        -3  3  0 ;
        -4  4  0 ;
        -5  5  0 ;
         0 -1  1]; % state change vectors dim jxn
propensity=@(x) [c0(1)*x(1)*x(2), c0(2)*x(1)*(x(1)-1)*x(2), c0(3)*x(1)*(x(1)-1)*(x(1)-2)*x(2),...
                 c0(4)*x(1)*(x(1)-1)*(x(1)-2)*(x(1)-3)*x(2), c0(5)*x(1)*(x(1)-1)*(x(1)-2)*(x(1)-3)*(x(1)-4)*x(2), c0(6)*x(2)]; 

simulation=SSA_mpd_2StopConditionALL(propensity,nu,x0,Tgrid,nsimula,Imax);
end
