function K=r3_matrix(m)
% Assembly only, before the stationary solver step; never solve target here.
s=mphmatrix(m,'sol1','out',{'K'},'rowscale','off',...
 'initmethod','sol','initsol','zero','extractafter','v1');
K=s.K;
end
