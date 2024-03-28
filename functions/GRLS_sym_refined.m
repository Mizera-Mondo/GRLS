function R = GRLS_sym_refined(f, alpha, beta, la_1)
%GRLS For order-1 Filter
[N, K] = size(f);
I = eye(N);
R = I;
Mdup = dupmat(N);
V1 = kron(f(:, 1)', I)*Mdup;
Q = zeros(size(V1'*V1));
q = zeros(size(V1'*f(:, 1)));

maxIter = 100;
maxLSIter = 100;
for k = 2:K
    
    x_ik = f(:, k);
    Vk = kron(f(:, k - 1)', I)*Mdup;
    Q = beta*Q + Vk'*Vk;
    q = beta*q + Vk'*x_ik;

    r = sym2vech(R);
    isConverge = false;
    isMaxIterReach = false;
    iter = 1;

    while ~isConverge && ~isMaxIterReach

        isArmijoReach = false;
        isLSMaxIterReach = false;
        targetFun = @(r) 0.5*r'*Q*r - r'*q + la_1*norm(r, 1);
        grad = Q*r - q + la_1*sign(r);
        ratio = alpha;
        iterLS = 1;

        while ~isArmijoReach && ~isLSMaxIterReach

            deltaR = -1*ratio*grad;
            expDescent = grad'*deltaR;
            realDescent = targetFun(r + deltaR) - targetFun(r);
            ratio = ratio/2;
            iterLS = iterLS + 1;
            isArmijoReach= realDescent < 0.3*expDescent;
            isLSMaxIterReach = iterLS > maxLSIter;

        end

        r = r + deltaR;
        iter = iter + 1;
        isConverge = norm(deltaR)/norm(r) < 1e-3;
        isMaxIterReach = iter > maxIter;

    end
       R = vec2squ(Mdup*r);
       
    if rand > 0.95
        disp(k);
    end
    
end
end
% cofShrink = 0.95;
% baseRatio = 10;
%     for tau = 1:tau_max
% %         for i = 1:M
% %         end
%         rh = sym2vech(R);
%         grad = Q*rh - q + la_1*sign(rh); % Line Search!!!
%         armijoReach = false;
%         maxIterReach = false;
%         lineSearchCount = 0;
%        
%         while ~armijoReach && ~maxIterReach
%             deltaRh = baseRatio*cofShrink^lineSearchCount*grad;
%             realDescent = targetFun(x_ik, Vk, rh - deltaRh) - targetFun(x_ik, Vk, rh);
%             expDescent = grad'*(-deltaRh);
%             if  realDescent < lineSearchThreshold*expDescent 
%                 armijoReach = true;
%             end
%             lineSearchCount = lineSearchCount + 1;
%             maxIterReach = lineSearchCount >= 100;
%         end
%         rh = rh - deltaRh;
%         R = vec2squ(Mdup*rh);
%         
%         



