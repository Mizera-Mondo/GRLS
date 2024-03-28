function R = GRLS_sym(f, alpha, beta, la_1, lineSearchThreshold, tau_max, t_max)
%GRLS For order-1 Filter
[N, K] = size(f);
I = eye(N);
R = rand(N, N);
Mdup = dupmat(N);
V1 = kron(f(:, 1)', I)*Mdup;
Q = V1'*V1;
q = V1'*f(:, 1);
cofShrink = 0.95;
baseRatio = 10;
targetFun = @(x_ik, Vk, rh) norm(x_ik - Vk*rh);
for k = 2:K
    
    x_ik = f(:, k);
    Vk = kron(f(:, k - 1)', I)*Mdup;
    Q = beta*Q + Vk'*Vk;
    q = beta*q + Vk'*x_ik;

    
    for tau = 1:tau_max
%         for i = 1:M
%         end
        rh = sym2vech(R);
        grad = Q*rh - q + la_1*sign(rh); % Line Search!!!
        armijoReach = false;
        maxIterReach = false;
        lineSearchCount = 0;
       
        while ~armijoReach && ~maxIterReach
            deltaRh = baseRatio*cofShrink^lineSearchCount*grad;
            realDescent = targetFun(x_ik, Vk, rh - deltaRh) - targetFun(x_ik, Vk, rh);
            expDescent = grad'*(-deltaRh);
            if  realDescent < lineSearchThreshold*expDescent 
                armijoReach = true;
            end
            lineSearchCount = lineSearchCount + 1;
            maxIterReach = lineSearchCount >= 100;
        end
        rh = rh - deltaRh;
        R = vec2squ(Mdup*rh);
        
        

    end
    if rand > 0.95
        disp(k);
    end
    
end
end
