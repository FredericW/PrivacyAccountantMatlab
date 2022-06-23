shift = 1;
shift_steps = ceil(n*shift);

subsampling = 0.5;
ps = [p ; zeros(shift_steps,1)];
qs = subsampling*[zeros(shift_steps,1); p]+(1-subsampling)*ps;


k = 10; %number of compositions

close all

logmgf = @(t) log(sum(ps.*((ps./(qs)+1e-30)).^t));
mgf = @(t) (sum(ps.*((ps./(qs)+1e-30)).^t));


% chebfun it
f = chebfun(@(t) k*logmgf(t)-log(t)-log(1+t), [1e-2,500],'splitting','on');
f_diff = diff(f);
f_hess = diff(f,2);

% "exact"
a = 0.01;
lims =  [-1500,1500];
f2 = chebfun(@(t) (mgf(a+1i*t)), lims,'splitting','on');

filter = chebfun(@(t) 1/((a+1i*t).*(a+1i*t+1)),lims,'splitting','on' );

t = chebfun('t',lims,'splitting','on');

% moments accountant

l_max = 100; %maximum lambda in MGF

fma = chebfun(@(t) logmgf(t), [0,l_max],'splitting','on'); %chebfun it!
l = chebfun(@(t) t, [0,l_max],'splitting','on'); % dummy variable

hold on
plot(f);
plot(f_diff);
%plot(f_hess);
grid on
set(gca, 'YScale', 'log')
hold off 
figure;

epsilon = transpose(0:0.1:5);
delta = zeros(length(epsilon),1);
delta_e = zeros(length(epsilon),1);
delta_m = zeros(length(epsilon),1);
% saddlepoint
for i=1:length(epsilon)
    length(epsilon)-i
    e = epsilon(i);
    x0 = roots(f_diff-e);

    delta(i) = exp(f(x0)-e*x0)/sqrt(abs(f_hess(x0)))/sqrt(2*pi);
    
    delta_e(i) = real(exp(-e*a)*sum((f2.^k)*exp(-e*(1i*t))*filter)/(2*pi));
    
    delta_m(i) = exp(min(k*fma-e*l));

end

hold on
plot(epsilon,delta,'b')
plot(epsilon,real(delta_e),'--r')
plot(epsilon,delta_m,'k')
set(gca, 'YScale', 'log')
xlabel('\epsilon')
ylabel('\delta')
legend({'Saddle point approximation','Numerical Integration','Moments Accountant'})
title(['Cactus accounting for ',num2str(k),'  compositions and subsampling rate ',num2str(subsampling)]);
set(gca,'FontSize',12)
grid on


hold off

saveas(gcf,'Saddle.png')
