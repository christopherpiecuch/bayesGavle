%%%%%%
% function model1: bayesian model for linear/parametric process, 
% perfect data, and without limiting or bounding dates
%%%%%%
function [ALPHA,BETA,GAMMA2]=model1

%%%%%%
% Say hello
pause(0.1), disp('Hello.  Things have started.'), pause(0.1)

%%%%%%
% Load data
load 20230628_data1.mat

%%%%%%%%%
% Set the seeds of the random number generator
rng(sum(clock))
%rng(234) % if reproducibility is needed

% define some parameters
K=numel(s); 

%%% number of draws to perform
NN_burn=1000; % warm-up draws
NN_post=10000; % post-warm-up draws
thin_period=10; % thin chains keeping 1 of 10
NN_burn_thin=NN_burn/thin_period; % Total number of burn-in to keep
NN_post_thin=NN_post/thin_period; % Total number of post-burn-in to keep
NN=NN_burn+NN_post; % Total number of draws to take 
NN_thin=NN_burn_thin+NN_post_thin;%  Total number of draws to keep

%%%%%%%%%
% Set hyperparameters
% specify uninformative priors
pp=[]; pp=polyfit(s,z,1);
mu=0; % slope prior mean:  0 is agnostic choice
zeta2=(10*abs(pp(1)))^2; % slope prior variance 
eta=0; % intercept prior mean: 0 is agnostic choice
sigma2=(10*abs(pp(2)))^2; % intercept prior variance
xi=3; % 
chi=1/3*var(z); % 
% note on xi and chi; see Tingley and Huybers (2010)
% https://doi.org/10.1175/2009JCLI3015.1 page 2778 their bullet point on
% the right column on "sigma^2" for an interpretation on the xi and chi
% values in terms of the shape of the prior inverse-gamma distribution
% placed on the variance parameter

%%%%%%%%%
% Allocate space for the sample arrays
ALPHA=zeros(NN,1);
BETA=zeros(NN,1);
GAMMA2=zeros(NN,1);

%%%%%%%%%
% Draw initial values
alpha=[]; alpha=0;
beta=[]; beta=0;
gamma2=[]; gamma2=min([1 1/randraw('gamma', [0,1/(chi),(xi)], [1,1])]);

%%%%%%%%%
% Loop through the Gibbs sampler
for nn=1:NN
    if mod(nn,100)==0
        disp([num2str(nn),' of ',num2str(NN),' iterations done.']), 
    end
    nn_thin=[]; nn_thin=ceil(nn/thin_period);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(alpha|.)
    V=mu/zeta2; PSI=1/zeta2;
    for k=1:K
        V=V+s(k)*(z(k)-beta)/gamma2;
        PSI=PSI+(s(k)^2)/gamma2;        
    end
    PSI=1/PSI;
    alpha=PSI*V+sqrt(PSI)*randn(1);
    clear V PSI
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(beta|.)
    V=eta/sigma2; PSI=1/(1/sigma2+K/gamma2);
    for k=1:K
        V=V+(z(k)-alpha*s(k))/gamma2;
    end
    beta=PSI*V+sqrt(PSI)*randn(1);
    clear V PSI

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(gamma2|.)
    SUM_K=0;
    for k=1:K
        DYKK=[]; DYKK=z(k)-alpha*s(k)-beta;
        SUM_K=SUM_K+DYKK^2;           
    end
   	gamma2=1/randraw('gamma', [0,1/(chi+1/2*SUM_K),...
     	(xi+K/2)], [1,1]);
   	clear SUM_K DYKK

    % update arrays
    ALPHA(nn)=alpha;
    BETA(nn)=beta;
    GAMMA2(nn)=gamma2;
end

%%%%%%
% Clear burnin
ALPHA(1:NN_burn)=[];
BETA(1:NN_burn)=[];
GAMMA2(1:NN_burn)=[];

%%%%%%
% Thin chains
ALPHA=ALPHA(thin_period:thin_period:NN_post);
BETA=BETA(thin_period:thin_period:NN_post);
GAMMA2=GAMMA2(thin_period:thin_period:NN_post);

return