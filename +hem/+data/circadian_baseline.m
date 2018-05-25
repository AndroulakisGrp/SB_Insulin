function [data_t, data_y] = circadian_baseline(var)
%HEM.DATA.CIRCADIAN_BASELINE Data from Jeremy's circadian paper
%   [DATA_T,DATA_Y] = HEM.DATA.CIRCADIAN_BASELINE(VAR) returns two vectors
%   representing the experimental data for the variable VAR.
%
%   Reference:
%   4.  Scheff JD, Calvano SE, Lowry SF, Androulakis IP: Modeling the
%       influence of circadian rhythms on the acute inflammatory response.
%       J Theor Biol 2010, 264(3):1068-1076.
%
%   See also HEM.DATA.HRV_LPS, HEM.PLOT.VARS

switch upper(var)
    case 'F'
        % Data: "Endogenous cortisol determines the circadian rhythm of
        % lipopolysaccharide- but not lipoteichoic acid-inducible
        % cytokine release" by Hermann, et al.
        data_t = [1, 5, 9, 13, 17, 21]; % hr
        data_y = [45, 70, 120, 70, 55, 50]; % ng/mL
        data_y = data_y/mean(data_y);
    case 'M'
        % Data: "Melatonin the 'light of night' in human biology and
        % adolescent idiopathic scoliosis" by Grivas and Savvidou
        data_t = [0, 1, 2, 8, 13, 18, 23]; % hr
        data_y = [29, 43, 40, 6, 4, 2, 9]; % pg/mL
        data_y = data_y/mean(data_y);
    case 'P'
        % Data: "Diurnal Rhythms of Pro-inflammatory Cytokines" by
        % Petrovsky, et al.
        % data_t = 0:23;
        % data_y = [130, 125, 132, 144, 166, 84, 80, 71, 41, 73, 81, 80, 92, 104, 66, 73, 76, 74, 91, 125, 140, 128, 111, 116]; % Relative levels (because it is LPS stimulated)
        % Outliers removed
        data_t = [0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23];
        data_y = [130, 125, 132, 144, 84, 80, 71, 73, 81, 80, 92, 104, 66, 73, 76, 74, 91, 125, 140, 128, 111, 116]; % Relative levels (because it is LPS stimulated)
        data_y = data_y/min(data_y);
    case 'A'
        % Data: "Diurnal Rhythmicity of Human Cytokine Production" by
        % Petrovsky, et al.
        data_t = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]; % hr
        data_y = [195, 146, 156, 184, 180, 91, 66, 67, 53, 57, 50, 51, 78, 78, 58, 53, 74, 63, 96, 66, 63, 120, 157, 224]; % Relative levels (because it is LPS stimulated)
        data_y = data_y/mean(data_y);
    case 'EPI'
        % Data: "Circadian immune measures in healthy volunteers" by
        % Kronfol, et al.
        data_t = 0:2:22; % hr
        data_y = [6, 12, 16, 9, 22, 24, 21, 19, 14, 16, 12, 11]; % pg/mL
        i = 4;
        j = 12;
        data_y = data_y/mean(data_y);
    case 'HRV'
        % Data: "Circadian rhythm of heart rate and heart rate variability"
        % by Massin et al.
        data_t = 1:24; % hr
        data_y = [152, 167, 184, 171, 194, 206, 168, 144, 95, 90, 122, 101, 80, 98, 97, 82, 130, 128, 123, 174, 168, 182, 195, 172]; % SDNN (ms)
        data_y = data_y/mean(data_y);
    otherwise
        data_t = [];
        data_y = [];
end
