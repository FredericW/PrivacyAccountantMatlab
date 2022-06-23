% They use 1e-3 subsampling, sigma = 0.8, and delta = 10^{-7}
%clear all
close all
% check if the is the correct PLRV
n = 250; % number of compositions
qs = 1e-3; % subsampling probability
sigma = .8; % variance of underlying Gaussian distribution
shift = 1; %shift

% these are the two PDFs for Poisson sub-sampling
gauss = @(t) 1/sqrt(2*pi)/sigma*exp(-t.^2/(2*sigma^2));
pq_ratio = @(t) 1 / ( (1-qs)+qs*(exp(( -(t-shift)^2+t^2  ) / (2*sigma^2))  )  ) ;

% range where pdfs are defined
range = [-11*sigma,11*sigma];

% create p distribution and PLRV el
p = chebfun(gauss,range,'splitting','on');
eL = chebfun(pq_ratio,range,'splitting','on');


% instantiate CGF and MGF
K_range = [1e-5,5000];

K = chebfun(@(t) log(sum(p*(eL)^t)),K_range,'splitting','on' ); % cgf
MGF = chebfun(@(t) (sum(p*(eL)^t)),K_range,'splitting','on' ); % mgf



% exact computation
as = 0.1; % shift
lims =  [-8000,8000];
te = chebfun('t',lims,'splitting','on');

f2 = chebfun(@(t) (sum(p*(eL)^(as+1i*t))), lims,'splitting','on');
filter = chebfun(@(t) 1/((as+1i*t).*(as+1i*t+1)),lims,'splitting','on' );

% functions for SPA
barrier = chebfun(@(t) -(1/t)-(1/(1+t)),K_range,'splitting','on');

d1K = diff(K); % first derivative
d2K = diff(K,2); % second derivative

mu = d1K(0); % mean (i.e., KL)

% compute third moment 
logeL = log(eL);

P = @(t) sum((p*abs(logeL-d1K(t))^3 )*(eL)^t)/MGF(t);


% q function
q = @(z) qfunc(z)*sqrt(2*pi)*z*exp(.5*z^2);

% gamma function
gamma = @(t,n,e) (n*d1K(t)-e)/sqrt(n*d2K(t));

alpha = @(t,n,g) sqrt(n*d2K(t))*t - g;

beta = @(t,n,g) sqrt(n*d2K(t))*(t+1) - g;

Delta = @(a,b,g) (q(a)/a - q(b)/b )*exp(-.5*g^2)/sqrt(2*pi);


% chebfun the root-finding function for SPA
root_fun = n*d1K + barrier;

% Edgeworth approximation
F = chebfun(@(t) n*K(t)-log(t)-log(1+t) ,K_range ,'splitting','on');
dF2 = diff(F,2);
dF3 = diff(F,3);
dF4 = diff(F,4);

% initialize vectors for plots/value tracking
max_epsilon = 5;
epsilon = transpose(0:0.005:max_epsilon);
delta_sp = zeros(length(epsilon),1);
error_sp = zeros(length(epsilon),1);
delta_e = zeros(length(epsilon),1);
delta_edge = zeros(length(epsilon),1);

% main iteration
for i=1:length(epsilon)
    length(epsilon)-i
    e = epsilon(i);

    % find saddle point
    t = roots(root_fun - e);
    
    % break if no saddle point found (need to increase chebfun lim)
    if size(t,1) == 0
        i = i-1;
        break
    end
   
    % compute quantities for SPA bound
    g = gamma(t,n,e);
    a = alpha(t,n,g);
    b = beta(t,n,g);

    d = Delta(a,b,g);

    ma = exp(n*K(t)-e*t);

    delta_sp(i) = ma*d;

    % error in SPA bound
    error_sp(i) = ma*exp(t*log(t)-(1+t)*log(1+t))*1.12*n*P(t)/( ( (n*d2K(t))^1.5 ) );

    % exact computation
    delta_e(i) = real(exp(-e*as)*sum((f2.^n)*exp(-e*(1i*te))*filter)/(2*pi));

    % edgeworth approximation
    approx_edge = exp(F(t)-e*t)/sqrt(2*pi*dF2(t)) ;
    error_edge = 1+ dF4(t)/ (8*dF2(t)^2 ) - (5*dF3(t)^2 )/(24*dF2(t)^3);
    delta_edge(i) = approx_edge*error_edge;

    % only plot for delta up to 1e-9
    if delta_e(i) <= 1e-9
        break
    end

end

epsilon = epsilon(1:i);
delta_e = delta_e(1:i);
delta_edge = delta_edge(1:i);
delta_sp = delta_sp(1:i);
error_sp = error_sp(1:i);

hold all
plot(epsilon,delta_e,'g')
plot(epsilon,delta_sp,'b')
plot(epsilon,delta_sp+error_sp,'k--')
plot(epsilon,delta_sp-error_sp,'k--')

plot(epsilon,delta_edge,'r')

set(gca, 'YScale', 'log')
xlabel('\epsilon')
ylabel('\delta')
legend({'Numerical integration','Saddle point approximation','SPA Upper Bound','SPA Lower Bound','Edgeworth'});
title(['Subsampled Gaussian accounting for ',num2str(n),'  compositions and subsampling rate ',num2str(qs)]);
set(gca,'FontSize',12)
grid on

saveas(gcf,'Saddle.png')

figure;

hold on
rel_error_sp = abs(100*(1-delta_sp./delta_e));
rel_error_edge = abs(100*(1-delta_edge./delta_e));

plot(epsilon,rel_error_sp,'b');
plot(epsilon,rel_error_edge,'r');

set(gca, 'YScale', 'log');
xlabel('\epsilon');
ylabel('\delta relative error (%)');
legend({'Saddle point approximation','Edgeworth'});
title(['Relative error to numerical integration for Subsampled Gaussian accounting'],['n=',num2str(n),' compositions and subsampling rate q=',num2str(qs)]);
set(gca,'FontSize',12)

saveas(gcf,'Relative error.png')
hold off


