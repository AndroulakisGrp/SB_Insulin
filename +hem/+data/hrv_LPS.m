function [data_t, data_y] = hrv_LPS(var)
%HEM.DATA.HRV_LPS Data from Pegy's HRV paper
%   [DATA_T,DATA_Y] = HEM.DATA.HRV_LPS(VAR) returns two vectors
%   representing the experimental data for the variable VAR.
%
%   Reference:
%   3.  Foteinou PT, Calvano SE, Lowry SF, Androulakis IP: Multiscale model
%       for the assessment of autonomic dysfunction in human endotoxemia.
%       Physiol Genomics 2010, 42(1):5-19.
%
%   See also HEM.DATA.CIRCADIAN_BASELINE, HEM.PLOT.VARS

switch upper(var)
    case 'F'
        data_t = [0, 1, 1.5, 2, 3, 4, 6, 24];
        data_y = [1.0000    1.1192    1.6606    2.0779    2.1378    2.0865    1.3548    1.1221];
    case 'P'
        data_t = [0, 2, 4, 6, 9, 24];
        data_y = [1.0000    7.9400    5.1600    2.7000    1.8700    1.1600];
    case 'A'
        data_t = [0, 2, 4, 6, 9, 24];
        data_y = [1.0000    1.2300    1.8800    2.1300    1.7000    0.9600];
    case 'EPI'
        data_t = [0, 2, 4 ,6];
        data_y = [1.0000    2.4842    1.6000    1.0211];
    case 'HRV'
        data_t = [0, 1, 2, 3, 4, 5, 6, 24];
        data_y = [1.0000    0.7826    0.4644    0.2975    0.3145    0.6043    0.5034    0.6859];
    otherwise
        data_t = [];
        data_y = [];
end
