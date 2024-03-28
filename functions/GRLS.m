function R = GRLS(f, alpha, beta, la_1)
%GRLS For order-1 Filter
% Ref: Graph Recursive Least Squares Filter for Topology Inference in Causal Data Processes
[N, K] = size(f);
I = eye(N);
R = I;
V1 = kron(f(:, 1)', I);
x_ik = f(:, 2);
Q = V1'*V1;
q = V1'*x_ik;

maxIter = 100;
maxLSIter = 100;
for k = 2:K - 1

    r = mat2vec(R);
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
    

%     if rand > 0.95
%         disp(k);
%     end

    R = vec2squ(r);
    x_ik = f(:, k + 1);
    Vk = kron(f(:, k)', I);
    Q = beta*Q + Vk'*Vk;
    q = beta*q + Vk'*x_ik;
end  
end
