%% GA
tic
clear; clc
load dataset
% 输入参数
% candidate_rs： 备选点坐标
% d：需求量
% best_restime：最佳救援时间 coverage
% ett：effective travel time distance_matrix
% 输出参数
% rs：被选择的站点编号

% 遗传算法参数设置
maxGenerations = 200; % 最大迭代次数
crossoverRate = 0.6; % 交叉率
mutationRate = 0.03; % 变异率

% 初始化结果存储
bestCoverageRates = zeros(length(candidate_rs), 1);
populationSize = 100; % 种群大小

% 遍历1到50个供应点个数
for numSupplyPoints = 1: length(candidate_rs)
    population = zeros(populationSize, length(candidate_rs));
    % 剩余的种群随机初始化
    for i = populationSize
        onesIdx = randperm(length(candidate_rs), numSupplyPoints);
        population(i, onesIdx) = 1;
    end
    
    
    bestCoverageRate = 0;
    fitnessHistory = zeros(maxGenerations, 1);
    coverageHistory = zeros(maxGenerations, 1);

    for generation = 1:maxGenerations
        % 计算适应度
        coverageRates = zeros(populationSize, 1);
        for i = 1:populationSize
            coverageRates(i) = calculateCoverageRate(population(i, :), ett, best_restime);
        end

        % 找到最佳覆盖率
        [maxCoverageRate, maxIdx] = max(coverageRates);
        if maxCoverageRate > bestCoverageRate
            bestCoverageRate = maxCoverageRate;
        end
        individual = population(maxIdx, :);
        rs = find(individual == 1);

        % 保存每代的最佳覆盖率值
        fitnessHistory(generation) = bestCoverageRate;
        coverageHistory(generation) = maxCoverageRate;

% %% 轮盘赌
% % 计算适应度总和
% totalFitness = sum(coverageRates);
% 
% % 检查是否存在适应度为负或零的情况，避免除以零或得到不合理的概率
% if totalFitness <= 0
%     error('Total fitness should be greater than zero.');
% end
% 
% % 计算每个个体的累积概率
% cumulativeProbabilities = coverageRates / totalFitness;
% cumulativeProbabilities = cumsum(cumulativeProbabilities);  % 使用 cumsum 函数计算累积概率
% 
% % 轮盘赌选择
% selectedIndices = zeros(1, populationSize);
% for i = 1:populationSize
%     % 生成一个[0,1]之间的随机数
%     randNum = rand();
%     % 找到随机数落在哪个区间
%     for j = 1:length(cumulativeProbabilities)
%         if randNum <= cumulativeProbabilities(j)
%             selectedIndices(i) = j;
%             break;  % 找到区间后退出循环
%         end
%     end
% end
% 
%        % 使用选中的个体索引更新种群
%        population = population(selectedIndices, :);
%%
        % 交叉
        for i = 1:2:populationSize
            if rand < crossoverRate
                crossoverPoint = randi(length(candidate_rs),1);
                offspring1 = [population(i, 1:crossoverPoint), population(i+1, crossoverPoint+1:end)];
                offspring2 = [population(i+1, 1:crossoverPoint), population(i, crossoverPoint+1:end)];
                population(i, :) = offspring1;
                population(i+1, :) = offspring2;
            end
        end

        % 变异
        for i = 1:populationSize
            for j = 1: length(candidate_rs)
                if rand < mutationRate
                    population(i, j) = 1 - population(i, j);
                end
            end
            % 确保变异后仍有 numSupplyPoints 个1
            if sum(population(i, :)) ~= numSupplyPoints
                onesIdx = find(population(i, :) == 1);
                zerosIdx = find(population(i, :) == 0);
                if length(onesIdx) > numSupplyPoints
                    % 多余的1变成0
                    toZero = randperm(length(onesIdx), length(onesIdx) - numSupplyPoints);
                    population(i, onesIdx(toZero)) = 0;
                elseif length(onesIdx) < numSupplyPoints
                    % 不足的0变成1
                    toOne = randperm(length(zerosIdx), numSupplyPoints - length(onesIdx));
                    population(i, zerosIdx(toOne)) = 1;
                end
            end
        end

        % 显示进度
        disp(['Supply Points ', num2str(numSupplyPoints), ', Generation ', num2str(generation), ': Best Coverage Rate = ', num2str(bestCoverageRate)]);

        % 判断是否达到最大代数
        if generation == maxGenerations
            break
        end
    end
    if bestCoverageRate == 1
        bestCoverageRates(numSupplyPoints) = bestCoverageRate;
        break
    end
    % 保存每个供应点个数下的最优解
    bestCoverageRates(numSupplyPoints) = bestCoverageRate;
end

% 绘制图表
% figure;
% hold on;
% plot(1: length(candidate_rs), bestCoverageRates, 'r*-');
% xlabel('Number of Supply Points');
% ylabel('Coverage Rate');
% title('Coverage Rate vs. Number of Supply Points');
% legend('Best Coverage Rates');
toc

