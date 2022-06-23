k = 1e4; % number of compositions
q = 0.01; % subsampling probability
sigma = 1; % variance of underlying Gaussian distribution
%eps = 1;
t = linspace(-5,5,1000)';
eps = linspace(0,1,100);


% these are the two PDFs for Poisson sub-sampling
fy = @(t) 1/sqrt(2*pi)/sigma*exp(-t.^2/(2*sigma^2));
fx = @(t) q*fy(t-1)+(1-q)*fy(t);

% this is non-zero to avoid the pole
a = 1;

% phi is the characteristic function
phi = chebfun(@(s) sum(chebfun(@(t) fx(t).*exp((a+1i.*s).*log(fx(t)./fy(t))),[-20,20])).^k,[0,20]);

% w is the function we need to integrate to compute delta
w = chebfun(@(s) real(phi(s).*exp(-(a+1i*s)*eps)./(a+1i*s)./(a+1+1i*s)),[0,20]);

% now do the integral
delta = sum(w)/pi;

plot(delta)
