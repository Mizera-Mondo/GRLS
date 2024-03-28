L = 200; % Length of signal
Nv = 30; % Num of verteces
Ne = 60; % Num of edges

betaStart = 0.98;
betaNum = 3;
betaEnd =  1.02;
betaStep = (betaEnd - betaStart)/(betaNum - 1);

laStart = 1;
laNum = 3;
laEnd = 2;
laStep = (laEnd - laStart)/(laNum - 1);
performanceCount = zeros(7, betaNum, laNum);

for count = 1:1

    disp('===================================================================');
    disp(['||Iteration: ' num2str(count)]);
    

    A = rand_digraph(Nv, Ne, 0.3, 0.1); % Adjacency, graph topology
    mu = 0.1*zeros(1, Nv); % Expectation of stimulus
    R = randn(Nv) + ones(Nv);
    sigma = triu(R)'*triu(R); % Covariance of stimulus
    S = repmat(mu, L, 1) + randn(L, Nv)*R; % Stimulus
    S = S';
    X = zeros(Nv, L);
    X(:, 1) = S(:, 1);
    for i = 2:L
        X(:, i) = A*X(:, i - 1) + S(:, i);
    end
    
    performanceCountCurrent = zeros(7, betaNum, laNum);
    
    betaCount = 1;

    for beta = betaStart:betaStep:betaEnd
        
        laCount = 1;

        for la_1 = laStart:laStep:laEnd

        tic;
        Aest = GRLS(X, 10, beta, la_1);
        T = toc;

        [acc, rec, pre, fM] = classifierPerformance(A > 0, Aest > 2e-1);
        performanceCountCurrent(:, betaCount, laCount) = [beta, la_1, acc, rec, pre, fM, T]';
        laCount = laCount + 1;

        end
        
        disp(['||Iteration ' num2str(count) ', ' num2str((betaCount/betaNum)*100) '% Completed']);
        betaCount = betaCount + 1;
        

    end
    performanceCount = performanceCount + performanceCountCurrent;
    disp(['||Iteration ' num2str(count) ' completed, saving current data']);
    save('gridSearchFineGRLS.mat', "performanceCount");
    disp('||Data saved.');
end