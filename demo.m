A = rand_digraph(30, 60, 0.4, 0.1); % 边的权重不宜太小，否则算法精度极差
L = 200;
stmls = 0.2*randn(30, L);
sig = stmls;
for i = 2:L
    sig(:, i) = A*sig(:, i - 1) + stmls(:, i);
end
Aest = GRLS(sig, 10, 1, 0.01); 
close all; 
imagesc(A); 
colorbar;
figure; 
imagesc(Aest);
colorbar;
err =@(A) norm(sig(:, 2:L) - A*sig(:, 1:L-1) - stmls(:, 2:L))/norm(stmls(:, 2:L));
disp(err(A));
disp(err(Aest));
% res = zeros(5, 1);
% for alpha = 0.1:0.2:1
%     for beta = 0.1:0.2:1
%         for la_1 = 0.1:0.2:1
%             for la_2 = 0.1:0.2:1
%                 Aest = GRLS(sig, alpha, beta, la_1, la_2, 10, 3);
%                 res = [res, [alpha, beta, la_1, la_2, norm(Aest - A)/norm(A)]'];
%             end
%         end
%     end
% end

% imagesc(Aest);
% figure; imagesc(A);
% Aest = GRLS(sig, 0.1, 0.1, 0.1, 0.1, 10, 10);