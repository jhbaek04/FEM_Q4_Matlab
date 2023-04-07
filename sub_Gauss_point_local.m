function [ GI ]  =  sub_Gauss_point_local
% hard-coded for 2-point rule.

GI = struct;

xi_1d = [ -1 ; 1 ] / sqrt(3);
wt_1d = [ 1 ; 1 ];


% For 2D,
% --------
% | 4  3 |
% | 1  2 |
% --------
xi_2d = [   [xi_1d(1); xi_1d(2); xi_1d(2); xi_1d(1)]  , ...
            [xi_1d(1); xi_1d(1); xi_1d(2); xi_1d(2)]   ];
wt_2d = [ 1 ; 1 ; 1 ; 1];


GI.xi_1d = xi_1d;
GI.weight_1d = wt_1d;

GI.xi_2d = xi_2d;
GI.weight_2d = wt_2d;

end