%% An algorithm based on greedy
tic
clear; clc
load dataset

% 初始化最佳方案和最佳覆盖率
best_solutions = [];
best_coverages = [];

for num_selected_points = 1: size(candidate_rs, 1)
    % 初始化种群
    % 从设置1个点开始
    if num_selected_points == 1
        population = eye(size(candidate_rs, 1));
        coverage = zeros(size(candidate_rs, 1), 1);
        for i = 1: length(population)
            coverage(i) = calculateCoverageRate(population(i, :), ett, best_restime);
        end
        index = find(coverage == max(coverage), 1);
        best_individual = population(index, :);
        best_solutions = [best_solutions; best_individual];
        best_coverages = [best_coverages, max(coverage)];
        disp(['当选择', num2str(num_selected_points), '个点时，最大覆盖率为', num2str(max(coverage))]);
    else
        population = zeros(size(candidate_rs, 1) - num_selected_points + 1, size(candidate_rs, 1));
        coverage = zeros(size(candidate_rs, 1) - num_selected_points + 1, 1);
        for i = 1: size(population, 1)
            population(i, :) = best_individual;
        end
        Index = find(best_individual ~= 1);
        for i = 1: length(Index)
            population(i, Index(i)) = 1;
        end
        for i = 1: size(population, 1)
            coverage(i) = calculateCoverageRate(population(i, :), ett, best_restime);
        end
        index = find(coverage == max(coverage), 1);
        best_individual = population(index, :);
        best_solutions = [best_solutions; best_individual];
        best_coverages = [best_coverages, max(coverage)];
        disp(['当选择', num2str(num_selected_points), '个点时，最大覆盖率为', num2str(max(coverage))]);
        if max(coverage) == 1
            break
        end
    end
end

for i = 1: length(best_coverages)
    disp(['当选择', num2str(i), '个点时']);
    disp('最佳选择方案为');
    index = find(best_solutions(i, :) == 1);
    disp(index);
    disp('最大覆盖率为');
    disp(best_coverages(i));
end
% selected_point = find(best_solutions(end, :) == 1);
toc
