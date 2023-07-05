%%%%%%
% function model4: bayesian model for nonlinear/nonparametric process, 
% uncertain data, and with limiting or bounding dates
% input is the form of covariance kernel you desire
%   kernNum=1 --> squared exponential
%   kernNum=2 --> matern v=3/2;p=1
%   kernNum=3 --> matern v=5/2;p=2
%   kernNum=4 --> % rational quadratic
%%%%%%
function [Y,TEE,GAMMA,LAMBDA,PHI,YPOST,YRATE,tpost,s,z,zet,sig,ind]=model4(kernNum)

%%%%%%
% Say hello
pause(0.1), disp('Hello.  Things have started.'), pause(0.1)

%%%%%%
% Load data
load 20230628_data4.mat

%%%%%%%%%
% Set the seeds of the random number generator
rng(sum(clock))
%rng(234) % if reproducibility is needed
 
%%% number of draws to perform
nnBurn=1000; % warm-up draws
nnPost=10000; % post-warm-up draws
nnThin=10; % thin chains keeping 1 of 10

%%%%%%%%%%%%%%%%%%%%%%%%
% sort data by age
[~,jj]=sort(s);
z=z(jj); s=s(jj); sig=sig(jj); zet=zet(jj); ind=ind(jj);
Sig=diag(sig.^2);
Zet=diag(zet.^2); 
invSig=inv(Sig);
invZet=inv(Zet);
J=numel(z);
clear data_*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what follows below is a slightly complicated procedure
% to select initial values for sea level, marine limits, and 
% terrestrial limits

% base prior values on optimized hyperparameters
% from the global holocene study by Khan et al. (2015)
% Curr Clim Change Rep (2015) 1:247–262
% doi 10.1007/s40641-015-0029-z
load('20230628_khanetal2015_tables1.mat')
etaphi=mean(phimle); chi2phi=var(phimle);
etagamma=mean(gammamle); chi2gamma=var(gammamle);
etalambda=mean(lambdamle); chi2lambda=var(lambdamle);
clear *mle

% before creating first-guess sea level, use maximum likelihood estimation 
% to find optimal parameter values
kk=find(ind==0);
fun = @(x)log(det(Zet(kk,kk)+exp(x(2))*eye(numel(kk))+exp(x(1))*setBayesKernal(abs(s(kk)-s(kk)'),x(3),kernNum)))+z(kk)'/(Zet(kk,kk)+exp(x(2))*eye(numel(kk))+exp(x(1))*setBayesKernal(abs(s(kk)-s(kk)'),x(3),kernNum))*z(kk);
x0 = [etagamma,etalambda,etaphi];
y = fminsearch(fun,x0);
gammamle=y(1);
lambdamle=max([log(0.1^2) y(2)]);
phimle=min([-6.9078 y(3)]);

% draw initial values
maxHeight=8800; % height of everest in m
minDepth=-11000; % depth of mariana in m
t=s;

aa=find(ind==0);
bb=1:numel(t);
T1=abs(t(aa)-t(aa)');
T2=abs(t(bb)-t(aa)');
w1=Zet(aa,aa)+exp(lambdamle)*eye(numel(aa))+exp(gammamle)*setBayesKernal(T1,phimle,kernNum);
w2=exp(gammamle)*setBayesKernal(T2,phimle,kernNum);
y=w2/w1*z(aa);
ind0=ind;
true=1;
while true
    ii=find((y>z&ind==1)|(y<z&ind==-1));
    if isempty(ii)
        true=0;
    end
    ind(ii)=0;
    aa=find(ind==0);
    bb=1:numel(t);
    T1=abs(t(aa)-t(aa)');
    T2=abs(t(bb)-t(aa)');
    w1=Zet(aa,aa)+exp(lambdamle)*eye(numel(aa))+exp(gammamle)*setBayesKernal(T1,phimle,kernNum);
    w2=exp(gammamle)*setBayesKernal(T2,phimle,kernNum);
    y=w2/w1*z(aa);
end
m=minDepth*ones(size(y));
f=maxHeight*ones(size(y));

clear aa bb T1 T2 w1 w2
ind=ind0; clear ind0
for jj=1:J
    if ind(jj)==-1
        m(jj)=z(jj)-1*zet(jj);
        if y(jj)<m(jj)
            y(jj)=m(jj)+2*zet(jj);
        end
    end
    if ind(jj)==1
        f(jj)=z(jj)+1*zet(jj);
        if y(jj)>f(jj)
            y(jj)=f(jj)-2*zet(jj);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% initial parameter values
phi=phimle;
gamma=[];
gamma=gammamle;
lambda=lambdamle;

%%%%%%%%%%%%%%%%%%%%%%%%
% other parameters 
N=nnBurn+nnPost;
ONE_J=ones(J,1);
ZERO_J=zeros(J,1);

%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output
Y=zeros(N,J);
M=zeros(N,J);
F=zeros(N,J);
TEE=zeros(N,J);
PHI=zeros(N,1);
GAMMA=zeros(N,1);
LAMBDA=zeros(N,1);

for nn=1:N
    if mod(nn,100)==0
        disp([num2str(nn),' of ',num2str(N),' iterations done.']), 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % define any matrices 
    T=abs(t-t');
    GMat=[]; GMat=exp(gamma)*setBayesKernal(T,phi,kernNum);
    LMat=[]; LMat=exp(lambda)*eye(J);
    wMat=GMat+LMat;
    invwMat=inv(wMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % evaluate y, f, m
    for jj=1:J %disp(num2str(jj))

        % evaluate y
        W=[]; Ws=[]; Wss=[]; invW=[];
        W=wMat; W(jj,:)=[]; W(:,jj)=[];
        invW=inv(W);
        Ws=wMat(jj,:); Ws(:,jj)=[];
        Wss=wMat(jj,jj);
        ytilde=y; ytilde(jj)=[];
        
        vy=[]; psiy=[]; II=[];
        II=ind(jj)==0;
        vy=II*z(jj)/(zet(jj))^2+inv(Wss-Ws*invW*Ws')*Ws*invW*ytilde;
        psiy=inv(II/(zet(jj))^2+inv(Wss-Ws*invW*Ws'));

        if psiy>=0
            sample=psiy*vy+sqrt(psiy)*randn(1);
            if sample>m(jj)&&sample<f(jj)
                y(jj)=sample; % only update if in range
            end
        end
        clear W Ws Wss invW ytilde vy psiy dummy sample
        
            % evaluate m
            if ind(jj)==-1 % marine limiting
                sample=z(jj)+zet(jj)*randn(1);
                if sample<y(jj)
                    m(jj)=sample; % only update if in range
                end
            else % not a marine limit; min
                m(jj)=minDepth;
            end
            
            % evaluate f
            if ind(jj)==1 % marine limiting
                sample=z(jj)+zet(jj)*randn(1);
                if sample>y(jj)
                     f(jj)=sample; % only update if in range
                end
            else % not a marine limit; min
                f(jj)=maxHeight;
            end
    end

	%%%%%%%%%%%%%%%%%%%%%%%%%
    % evaluate t
    t_now=t;
    t_cov=0.1*Sig;
    t_prp=mvnrnd(t_now,t_cov)';
    T_now=abs(t_now-t_now');
    T_prp=abs(t_prp-t_prp');
    w_now=LMat+exp(gamma)*setBayesKernal(T_now,phi,kernNum);
    w_prp=LMat+exp(gamma)*setBayesKernal(T_prp,phi,kernNum);
    invw_now=inv(w_now);
    invw_prp=inv(w_prp);
    ins_now=-0.5*(s-t_now)'*invSig*(s-t_now)-0.5*y'*invw_now*y;
    ins_prp=-0.5*(s-t_prp)'*invSig*(s-t_prp)-0.5*y'*invw_prp*y;
    MetFrac=det(w_prp*invw_now)^(-1/2)*exp(ins_prp-ins_now);
    success_rate=min(1,MetFrac);
    if rand(1)<=success_rate
        t_now=t_prp;
    end
    t=t_now;
    T=abs(t-t');
    GMat=[]; GMat=exp(gamma)*setBayesKernal(T,phi,kernNum);
    LMat=[]; LMat=exp(lambda)*eye(J);
    wMat=GMat+LMat;
    invwMat=inv(wMat);
    
    if nn>((nnBurn)/2)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(phi|.)
    Phi_now=phi;
    Phi_std=2;
    Phi_prp=normrnd(Phi_now,Phi_std);
    w_now=LMat+exp(gamma)*setBayesKernal(T,Phi_now,kernNum);
    w_prp=LMat+exp(gamma)*setBayesKernal(T,Phi_prp,kernNum);
    invw_now=inv(w_now);
    invw_prp=inv(w_prp);
    ins_now=-1/(2*chi2phi)*(Phi_now-etaphi)^2-1/2*y'*invw_now*y;
    ins_prp=-1/(2*chi2phi)*(Phi_prp-etaphi)^2-1/2*y'*invw_prp*y;
    MetFrac=det(w_prp*invw_now)^(-1/2)*exp(ins_prp-ins_now);
    success_rate=min(1,MetFrac);
    if rand(1)<=success_rate
        Phi_now=Phi_prp;
        success(nn)=1;
    else
        success(nn)=0;
    end
    phi=Phi_now;
    T=abs(t-t');
    GMat=[]; GMat=exp(gamma)*setBayesKernal(T,phi,kernNum);
    LMat=[]; LMat=exp(lambda)*eye(J);
    wMat=GMat+LMat;
    invwMat=inv(wMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(gamma|.)
    Gam_now=gamma;
    Gam_std=2;
    Gam_prp=normrnd(Gam_now,Gam_std);
    w_now=LMat+exp(Gam_now)*setBayesKernal(T,phi,kernNum);
    w_prp=LMat+exp(Gam_prp)*setBayesKernal(T,phi,kernNum);
    invw_now=inv(w_now);
    invw_prp=inv(w_prp);
    ins_now=-1/(2*chi2gamma)*(Gam_now-etagamma)^2-1/2*y'*invw_now*y;
    ins_prp=-1/(2*chi2gamma)*(Gam_prp-etagamma)^2-1/2*y'*invw_prp*y;
    MetFrac=det(w_prp*invw_now)^(-1/2)*exp(ins_prp-ins_now);
    success_rate=min(1,MetFrac);
    if rand(1)<=success_rate
        Gam_now=Gam_prp;
        success2(nn)=1;
    else
        success2(nn)=0;
    end
    gamma=max([log(.1^2) Gam_now]);
    T=abs(t-t');
    GMat=[]; GMat=exp(gamma)*setBayesKernal(T,phi,kernNum);
    LMat=[]; LMat=exp(lambda)*eye(J);
    wMat=GMat+LMat;
    invwMat=inv(wMat);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample from p(lambda|.)
    Lam_now=lambda;
    Lam_std=1;
    Lam_prp=normrnd(Lam_now,Lam_std);
    w_now=exp(Lam_now)*eye(J)+exp(gamma)*setBayesKernal(T,phi,kernNum);
    w_prp=exp(Lam_prp)*eye(J)+exp(gamma)*setBayesKernal(T,phi,kernNum);
    invw_now=inv(w_now);
    invw_prp=inv(w_prp);
    ins_now=-1/(2*chi2lambda)*(Lam_now-etalambda)^2-1/2*y'*invw_now*y;
    ins_prp=-1/(2*chi2lambda)*(Lam_prp-etalambda)^2-1/2*y'*invw_prp*y;
    MetFrac=det(w_prp*invw_now)^(-1/2)*exp(ins_prp-ins_now);
    success_rate=min(1,MetFrac);
    if rand(1)<=success_rate
        Lam_now=Lam_prp;
        success3(nn)=1;
    else
        success3(nn)=0;
    end
    lambda=max([log(0.01^2) Lam_now]);
    
    end % only start sampling parameters after half of burn in (viz. tingley and huybers 2010)

    % define output
    Y(nn,:)=y;
    M(nn,:)=m;
    F(nn,:)=f;
    TEE(nn,:)=t;
    GAMMA(nn)=gamma;
    PHI(nn)=phi;
    LAMBDA(nn)=lambda;

end

% delete burn in
Y(1:nnBurn,:)=[];
F(1:nnBurn,:)=[];
M(1:nnBurn,:)=[];
TEE(1:nnBurn,:)=[];
GAMMA(1:nnBurn)=[];
PHI(1:nnBurn)=[];
LAMBDA(1:nnBurn)=[];
N=N-nnBurn;

% thin chains
Y=Y(1:nnThin:end,:);
F=F(1:nnThin:end,:);
M=M(1:nnThin:end,:);
TEE=TEE(1:nnThin:end,:);
GAMMA=GAMMA(1:nnThin:end);
PHI=PHI(1:nnThin:end);
LAMBDA=LAMBDA(1:nnThin:end);
N=numel(LAMBDA);

% posterior prediction to get things on a regular time grid
tpost=0:100:1.2e4;
%tpost=0:100:tau;
KP=numel(tpost);
ONE_P=ones(KP,1);
YPOST=zeros(N,KP);
YRATE=zeros(N,KP);
for nn=1:N %disp(num2str(nn))
    
    %%%%%%%%%%%%%%%%%%%
    % smoothed process on regular intervals
    f=[]; f=Y(nn,:)';
    ttemp=[]; ttemp=abs(TEE(nn,:)'-TEE(nn,:));
    Kxx=exp(LAMBDA(nn))*eye(size(ttemp))+exp(GAMMA(nn))*setBayesKernal(ttemp,PHI(nn),kernNum);
    ttemp=[]; ttemp=abs(TEE(nn,:)'-tpost);
    Kxxs=exp(GAMMA(nn))*setBayesKernal(ttemp,PHI(nn),kernNum);
    ttemp=[]; ttemp=abs(tpost'-tpost);
    Kxsxs=exp(GAMMA(nn))*setBayesKernal(ttemp,PHI(nn),kernNum)+1e-4*eye(size(ttemp));

    mat=[]; mat=Kxsxs-Kxxs'/Kxx*Kxxs;
    mat=0.5*(mat+mat');
    
    fs=[]; fs=mvnrnd(Kxxs'/Kxx*f,mat)';
    YPOST(nn,:)=fs;  

	%%%%%%%%%%%%%%%%%%%
    % rate on regular intervals
    f=[]; f=Y(nn,:)';
    ttemp=[]; ttemp=abs(TEE(nn,:)'-TEE(nn,:));
    Kxx=exp(LAMBDA(nn))*eye(size(ttemp))+exp(GAMMA(nn))*setBayesKernal(ttemp,PHI(nn),kernNum);
    ttemp=[]; ttemp=(TEE(nn,:)'-tpost);
    Kxxs=exp(GAMMA(nn))*setBayesKernal01(ttemp,PHI(nn),kernNum);
    ttemp=[]; ttemp=abs(tpost'-tpost);
    Kxsxs=exp(GAMMA(nn))*setBayesKernal11(ttemp,PHI(nn),kernNum)+1e-10*eye(size(ttemp));
    mat=[]; mat=Kxsxs-Kxxs'/Kxx*Kxxs;
    mat=0.5*(mat+mat');
 
    fs=[]; fs=mvnrnd(Kxxs'/Kxx*f,mat)';
    YRATE(nn,:)=-fs;  
      
end