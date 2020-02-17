%% Load data

pr_fit_hr = dlmread("14K_fitness_highres.txt")';
pr_fit = dlmread("14K_fitness.txt")';
pr_param = dlmread("14K_params.txt");
mc1_fit = dlmread("14K_monteCarlo-1000_fitness.txt");
mc2_fit = dlmread("14K_monteCarlo-4000_fitness.txt");
mc1_param = dlmread("14K_monteCarlo-1000_params.txt");
mc2_param = dlmread("14K_monteCarlo-4000_params.txt");

%% Simple pareto front
figure;
scatter(pr_fit(:,2), pr_fit(:, 1), 25,'filled', 'MarkerEdgeColor', [0, 0, 0])
hold on
p3 = scatter([pr_fit(105, 2)], [pr_fit(105, 1)], 'xg');


grid on
legend("Pareto front", "Selected individual")
xlabel("Maximum GDOP [-]")
ylabel("Mean GDOP [-]")
% xlim([0, 1.5e4])
% title("Comparison of Pareto frontier with Monte-Carlo solutions")


%% Compare optimization results to monte-carlo to check convergence

figure;
scatter(pr_fit(:,2), pr_fit(:, 1), 25,'filled', 'MarkerEdgeColor', [0, 0, 0])
hold on
scatter([mc1_fit(2, :), mc2_fit(2, :)], [mc1_fit(1, :), mc2_fit(1, :)], '.', 'MarkerFaceColor',[.7 .07 .0])


legend("Pareto front", "Monte-Carlo")
grid on
xlabel("Maximum GDOP [-]")
ylabel("Mean GDOP [-]")
xlim([0, 1e5])
title("Comparison of Pareto frontier with Monte-Carlo solutions")

%% Compare optimization results to high resolution re-runs to check sensitivity to the time spacing.

figure;
hold on;
for i = 1:size(pr_fit, 1)
   plot([pr_fit(i, 2), pr_fit_hr(i, 2)], [pr_fit(i, 1), pr_fit_hr(i, 1)] , 'color', [0, 0, 0])
end
p1 = scatter(pr_fit(:,2), pr_fit(:, 1), 25,'filled', 'MarkerFaceColor',[0 0.4470 0.7410], 'MarkerEdgeColor', [0, 0, 0]);
p2 = scatter(pr_fit_hr(:,2), pr_fit_hr(:, 1), 25,'filled', 'MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor', [0, 0, 0]);
p3 = scatter([pr_fit(105, 2), pr_fit(60, 2)], [pr_fit(105,1), pr_fit(60, 1)], 'xg');

legend([p1 p2 p3], {'dt = 200', 'dt = 10', 'selected'})
grid on
xlabel("Maximum GDOP [-]")
ylabel("Mean GDOP [-]")
title("Sensitivity to solution resolution")

hold off;
magnifyOnFigure;
disp('Press a key...')
pause;

% Select number: 105 (actual solution, inclined); 60 (comparison, low_inc)


%% Plot design parameter plots

figure;
subplot(2, 2, 1)
scatter(rad2deg((pr_param(:,1)-0.5) * (pi/2)), pr_fit(:,1), 25,'filled', 'MarkerEdgeColor', [0, 0, 0])
hold on
scatter([rad2deg((pr_param(105,1)-0.5) * (pi/2)), rad2deg((pr_param(60,1)-0.5) * (pi/2))], [pr_fit(105,1), pr_fit(60, 1)], 'xg')
% scatter(rad2deg((mc1_param(1, 589)-0.5) * (pi/2)), mc1_fit(1, 589), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])
% scatter(rad2deg((mc2_param(1, 2144)-0.5) * (pi/2)), mc2_fit(1, 2144), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])
disp(pr_param(105,:))
disp(pr_param(60, :))


grid on
ylabel("Mean GDOP [-]")
xlabel("\Deltai_1 [deg]")

subplot(2, 2, 2)
scatter(rad2deg((pr_param(:,2)-0.5) * (pi/2)), pr_fit(:,1), 25,'filled', 'MarkerEdgeColor', [0, 0, 0])
hold on
scatter([rad2deg((pr_param(105,2)-0.5) * (pi/2)), rad2deg((pr_param(60,2)-0.5) * (pi/2))], [pr_fit(105,1), pr_fit(60, 1)], 'xg')
% scatter(rad2deg((mc1_param(2, 589)-0.5) * (pi/2)), mc1_fit(1, 589), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])
% scatter(rad2deg((mc2_param(2, 2144)-0.5) * (pi/2)), mc2_fit(1, 2144), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])

grid on
ylabel("Mean GDOP [-]")
xlabel("\Deltai_2 [deg]")

subplot(2, 2, 3)
scatter(rad2deg(pr_param(:,3)* (2/3) * pi), pr_fit(:,1), 25,'filled', 'MarkerEdgeColor', [0, 0, 0])
hold on
scatter([rad2deg(pr_param(105,3)* (2/3) * pi), rad2deg(pr_param(60,3)* (2/3) * pi)], [pr_fit(105,1), pr_fit(60, 1)], 'xg')
% scatter(rad2deg((mc1_param(3, 589)-0.5) * (pi/2)), mc1_fit(1, 589), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])
% scatter(rad2deg((mc2_param(3, 2144)-0.5) * (pi/2)), mc2_fit(1, 2144), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])

grid on
ylabel("Mean GDOP [-]")
xlabel("\Delta\theta_2 [deg]")

subplot(2, 2, 4)
scatter(rad2deg(pr_param(:,4) * (2/3) * pi), pr_fit(:,1), 25,'filled', 'MarkerEdgeColor', [0, 0, 0])
hold on
scatter([rad2deg(pr_param(105,4)* (2/3) * pi), rad2deg(pr_param(60,4)* (2/3) * pi)], [pr_fit(105,1), pr_fit(60, 1)], 'xg')
% scatter(rad2deg((mc1_param(4, 589)-0.5) * (pi/2)), mc1_fit(1, 589), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])
% scatter(rad2deg((mc2_param(4, 2144)-0.5) * (pi/2)), mc2_fit(1, 2144), 25,'filled', 'MarkerFaceColor',[.07 .7 .0], 'MarkerEdgeColor', [0, 0, 0])

grid on
ylabel("Mean GDOP [-]")
xlabel("\Delta\theta_3 [deg]")
