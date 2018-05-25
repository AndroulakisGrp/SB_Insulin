function [z, mus] = bifurcation(noise_to_compartment)

tic;

niter = 20; % Number of simulations at each noise level

% Noise
dmu = 0.05; % Step size
mus = [0.001, 0.01, .05:.1:2.15];

nmu = length(mus); % Number of noise levels

z = zeros(nmu, niter);

if matlabpool('size') == 0
    matlabpool open 8;
end

system('cat /dev/null > progress.txt'); % Clear progress file. LINUX ONLY!
parfor j=1:nmu
    s = struct;
    s.t_pre = 48;
    s.t_final = 24;
    s.stochastic_ensemble = true;
    s.mu = mus(j);
    s.circadian = true;
    s.noise_to_compartment = noise_to_compartment;
    qqq = zeros(1, niter);
    for i=1:niter
        [t, y] = hem.model.run(0.8, 9, s); % 0.9 is closer
        qqq(i) = y.P(end);
    end
%     y = run_model_normal(make_plots, niter, mus(j));

    z(j,:) = qqq;
    
    % Write a line to progress.txt and then count how many lines there are
    % to assess progress. LINUX ONLY!
    system('echo ''1'' >> progress.txt');
    [status, result] = system('cat progress.txt | wc -l');
    fprintf('%2.0f%%...\n', str2double(result)/nmu*100);
end
fprintf('\n');


%%
bins = 0:1:16;
% bifurcation2_heatmap(mus, z, niter, bins);
% %% Histogram plot
% figure;
% hold on;
% x = linspace(min(z(:))-1, max(z(:))+1, 40);
% x = linspace(min(z(:))-1, 25, 40);
% dmu_equiv = 0;
% for j=1:nmu
%     n = hist(z(j,:), x);
%     basex = mus(j)*ones(size(n));
    
    % How wide should the histograms be? Make the width equal to the
    % biggest peak in the smallest noise level. This might not always work.
%     if dmu_equiv == 0
%         dmu_equiv = dmu/max(n);
%     end
% %     xx = basex + n*dmu/niter;
%     xx = basex + n*dmu_equiv;
% %     plot(xx, x);
%     patch([mus(j), xx, mus(j), mus(j)], [x(1), x, x(end), x(1)], 'r', 'LineWidth', 1, 'EdgeColor', [1 0 0]);
%     plot(basex, x, 'k', 'LineWidth', 1);
% end
% xlim([0.04, mus(end)+dmu])
% ylim([min(x), max(x)]);
% xlabel('Noise level');
% ylabel('Steady state');
% set(gca, 'xtick', mus(1:3:end));

toc;

%%
figure;
bifurcation2_heatmap(mus, z, niter, bins, 1000, true);
% my_xticklabels_image([1, 2, 3:2:length(mus)], {'10^{-3}', '10^{-2}', '0.05', '0.15', '0.25', '0.35', '0.45'});
set(gca, 'XTick', 1:3:length(mus));
set(gca, 'XTickLabel', mus(1:3:length(mus)));
