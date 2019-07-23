function f = colebrook(Re,epsilon)
%COLEBROOK Computes the friction factor in pipes for given values of the 
%   Reynolds number (Re) and the relative roughness coefficient (epsilon).
%
%   Syntax:
%      f = colebrook(Re,epsilon)
%
%   Example 1: Single Re, single epsilon
%      Re = 1e5;
%      epsilon = 1e-4;
%      f = colebrook(Re,epsilon)
%
%   Example 2: Multiple Re, single epsilon
%      Re = 5000:1000:100000;
%      epsilon = 1e-4;
%      f = colebrook(Re,epsilon*ones(size(Re)));
%      plot(Re,f)
%
%   Example 3: Multiple Re, multiple epsilon
%      Re = logspace(4,8,100);
%      epsilon = linspace(1e-4,1e-1,100);
%      [RE,EPSILON] = meshgrid(Re,epsilon);
%      F = colebrook(RE,EPSILON);
%      surf(RE,EPSILON,F)
%
%   References:
%      [1] Colebrook, C. F., & White, C. M. (1937). Experiments with fluid
%          friction in roughened pipes. Proceedings of the Royal Society of
%          London. Series A - Mathematical and Physical Sciences, 161(906),
%          367-381.
%      [2] Colebrook, C. (1939). TURBULENT FLOW IN PIPES, WITH PARTICULAR
%          REFERENCE TO THE TRANSITION REGION BETWEEN THE SMOOTH AND ROUGH
%          PIPE LAWS. Journal of the Institution of Civil Engineers, 11(4),
%          133-156.
%
%   Author:
%      Ildeberto de los Santos Ruiz
%      idelossantos@ittg.edu.mx

shape = size(Re);
Re = Re(:);
epsilon = epsilon(:);
f = zeros(size(Re));
for k = 1:numel(Re)
    f(k) = fzero(@(f) 1/sqrt(f)+2*log10(epsilon(k)/3.7+...
        2.51/(Re(k)*sqrt(f))),[eps,1]);
end
f = reshape(f,shape);
end