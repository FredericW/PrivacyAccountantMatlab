clear all
close all

mu = [0;2;-2]; % means
sigma = [1]; % variance
gm = gmdistribution(mu,sigma); % mixture distribution

% plot mixture distribution
x = linspace(-5,5,1000)';
plot(x,pdf(gm,x));
title('Additive noise distribution');
ylabel('PDF');

d = 1; % shift (i.e., "neighboring" defn)

% uncomment to plot the function g that produces the log-loss r.v.
%y = g(x,d,gm);

% hold on
% plot(x,y);
CF=  chebfun (@(s) sum(chebfun (@ (x) pdf(gm,x).*exp(1i.*x.*s),[-20,20])), [0,20]);
figure;
s = linspace(0,10,1000);
plot(s,real(CF(s)));

% tic;
% z = random(gm,1e8); % draw many samples from the mixute
% toc;
% hold off
% figure;
% histogram(g(z,d,gm),2000,'DisplayStyle','stairs','Normalization','pdf'); % plot hist of transformed r.v.
% title('Privacy Loss Distribution (Empirical)');
% ylabel('Empirical PDF');

% function for computing privacy loss r.v.
function y = g(x,d,gm)
    y = log(pdf(gm,x)./pdf(gm,x+d));
end