function circadian()

tlpses = 0:0.5:23;

s.t_pre = 48;
s.paramset = 'nonnegative';

for i=1:length(tlpses)
    [t, y] = hem.model.run(.2, tlpses(i), s);
%     figure;
%     hem.plot.vars(t, y);
%     keyboard
    max_p(i) = max(y.P) % P
    max_epi(i) = max(y.EPI) % EPI
    min_hrv(i) = min(y.HRV) % HRV (min!!!)
end


%% PLOTTING
figure;
[AX, H1, H2] = plotyy(tlpses, min_hrv, tlpses, max_p);

set(H1, 'LineWidth', 2);
set(H2, 'LineWidth', 2);

set(get(AX(1), 'Ylabel'), 'String', 'Maximum depression in HRV (HRV_{min})');
set(get(AX(2), 'Ylabel'), 'String', 'Strength of the inflammatory response (P_{max})');
% set(get(AX(2), 'Ylabel'), 'String', 'EPI_{max}');

set(AX,'Box','on');
set(AX(1),'XTick', 0:6:24);
set(AX(1), 'XTickLabel', {'12am', '6am', '12pm', '6pm', '12am'});
set(AX(2), 'XTick', []);
xlim([0, 24]);
set(AX, 'xlim', [0, 24]);
grid on;
xlabel('Time of inflammatory stimulus (Hr)');
%ylabel('Strength of the inflammatory response (P_{max})');
%ylabel('Maximum depression in HRV (HRV_{min})');