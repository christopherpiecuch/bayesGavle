%
%
%
function K=setBayesKernal01(T,phi,kernalNumber)
if kernalNumber==1 % squared exponential
    %K=exp(-0.5*(exp(phi)*T).^2);
    K=(exp(phi)^2)*T.*exp(-0.5*(exp(phi)*T).^2);
elseif kernalNumber==2 % matern v=3/2;p=1
    %K=(ones(size(T))+sqrt(3)*(exp(phi)*T)).*exp(-sqrt(3)*(exp(phi)*T));
    K=(3*(exp(phi)^2*T)).*exp(-sqrt(3)*(exp(phi)*abs(T)));
elseif kernalNumber==3 % mater v=5/2;p=2
    %K=(ones(size(T))+sqrt(5)*(exp(phi)*T)+(5/3)*(exp(phi)*T).^2).*exp(-sqrt(5)*(exp(phi)*T));
    K=(5/3)*(exp(phi)^2)*(ones(size(T))+sqrt(5)*exp(phi)*abs(T)).*exp(-sqrt(5)*(exp(phi)*abs(T))).*T;
elseif kernalNumber==4 % rational quadratic
    alfa=1;
    %K=(ones(size(T))+1/(2*alfa)*(exp(phi)*T).^2).^(-alfa);
    K=(exp(phi)^2)*T.*(ones(size(T))+1/(2*alfa)*(exp(phi)*T).^2).^(-alfa-1);
end
return