% --------------------------------------------------------------------
% This MATLAB script plots two figures:
% (1) The figure on the left compares the block error rate (BLER) of
%     symbol-level GRAND with that of bit-level GRAND and uncoded 16-QAM.
% (2) The figure on the right compares the average number of error patterns
%     tested by the decoding algortihms until an error pattern was
%     accepted.
%
% The figure appears in:
% "Symbol-Level GRAND for High-Order Modulation over Block Fading Channels",
% IEEE Communications Letters, vol. 27, no. 2, Feb. 2023, pp. 447-451.
% 
% Authors: 
% Ioannis Chatzigeorgiou, Lancaster University, United Kingdom
% Francisco Monteiro, Instituto de Telecommunicacoes 
% and ISCTE-Instituto Universitario de Lisboa, Portugal
%
% If you use the simulation results of symbol-level GRAND, 
% please cite the paper above.
%
% ** Note that MATLAB for MacOS was used to obtain this figure. 
% The figure might not display properly on other platforms. **
% --------------------------------------------------------------------

% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
  
pos_figure1 = [0 0 500 260];
set(figure1,'Position',pos_figure1)

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.101785714285714 0.1699604743083 0.380357142857143 0.798418972332016]);

% Plot simulation results for uncoded 16-QAM (AWGN)
load('BLER_Rayleigh_uncoded_16QAM.mat', 'Eb_N0', 'qam16_BLER');
semilogy(Eb_N0, qam16_BLER, 'k-', 'linewidth', 1,'Parent',axes1);
hold on
clear qam16_BLER Eb_N0

% Create text
text('Parent',axes1,'FontSize',12,'Rotation',-37.5,'String','Uncoded 16-QAM',...
    'Position',[26.4210525387212 0.0226865827731463 0]);

% Plot simulation results for RLC + 16-QAM + Bit-level GRAND (Rayleigh, w_th=2)
load('BLER_Rayleigh_RLC_16QAM_bit_GRAND_mw2.mat', 'Eb_N0', 'qam16_BLER');
semilogy(Eb_N0(1:2:end), qam16_BLER(1:2:end), 'rs--', 'linewidth', 1, 'MarkerSize', 10, 'Parent', axes1);
hold on
clear qam16_BLER Eb_N0
% Plot simulation results for RLC + 16-QAM + Symbol-level GRAND (Rayleigh, w_th=2)
load('BLER_Rayleigh_RLC_16QAM_symbol_GRAND_mw2tr3.mat', 'Eb_N0', 'qam16_BLER');
semilogy(Eb_N0(1:2:end), qam16_BLER(1:2:end), 'bo--', 'linewidth', 1, 'Parent', axes1);
hold on
clear qam16_BLER Eb_N0
% Plot simulation results for RLC + 16-QAM + Bit-level GRAND (Rayleigh, w_th=3)
load('BLER_Rayleigh_RLC_16QAM_bit_GRAND_mw3.mat', 'Eb_N0', 'qam16_BLER');
semilogy(Eb_N0(2:2:end), qam16_BLER(2:2:end), 'rs-', 'linewidth', 1, 'MarkerSize', 10, 'Parent', axes1);
hold on
clear qam16_BLER Eb_N0
% Plot simulation results for RLC + 16-QAM + Symbol-level GRAND (Rayleigh, w_th=3)
load('BLER_Rayleigh_RLC_16QAM_symbol_GRAND_mw3tr5.mat', 'Eb_N0', 'qam16_BLER');
semilogy(Eb_N0(1:2:end), qam16_BLER(1:2:end), 'bo-', 'linewidth', 1, 'Parent', axes1);
hold on
clear qam16_BLER Eb_N0

xlim(axes1,[20 34]);
ylim(axes1,[10^(-3) 10^(-1)]);
set(axes1,'xtick',20:2:34, 'FontSize',14)
box(axes1,'on');
grid(axes1,'on');

ylabel('Block error rate (BLER)','FontSize',14, 'Parent', axes1);
xlabel('$E_b/N_0$ (dB)', 'Interpreter','latex','FontSize',16, 'Parent', axes1);
  
% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.59642857142857 0.1699604743083 0.380357142857143 0.798418972332016]);

semilogy(-1,1,'rs','MarkerSize', 10);
hold on
semilogy(-1,1,'bo');
hold on

% Plot simulation results for RLC + 16-QAM + Bit-level GRAND (w=2)
load('BLER_Rayleigh_RLC_16QAM_bit_GRAND_mw2.mat', 'Eb_N0', 'avg_patterns_tested');
semilogy(Eb_N0, avg_patterns_tested, 'rs--', 'linewidth', 1, 'MarkerSize', 10, 'Parent', axes2);
hold on;
%avg_patterns_tested_bit_GRAND = avg_patterns_tested;
clear avg_patterns_tested Eb_N0

% Plot simulation results for RLC + 16-QAM + Symbol-level GRAND (w=2)
load('BLER_Rayleigh_RLC_16QAM_symbol_GRAND_mw2tr3.mat', 'Eb_N0', 'avg_patterns_tested');
semilogy(Eb_N0, avg_patterns_tested, 'bo--', 'linewidth', 1, 'Parent', axes2);
hold on;
%avg_patterns_tested_symbol_GRAND1 = avg_patterns_tested;
clear avg_patterns_tested Eb_N0

% Plot simulation results for RLC + 16-QAM + Bit-level GRAND (w=3)
load('BLER_Rayleigh_RLC_16QAM_bit_GRAND_mw3.mat', 'Eb_N0', 'avg_patterns_tested');
semilogy(Eb_N0, avg_patterns_tested, 'rs-', 'linewidth', 1, 'MarkerSize', 10, 'Parent', axes2);
hold on;
%avg_patterns_tested_bit_GRAND = avg_patterns_tested;
clear avg_patterns_tested Eb_N0

% Plot simulation results for RLC + 16-QAM + Symbol-level GRAND (w=3)
load('BLER_Rayleigh_RLC_16QAM_symbol_GRAND_mw3tr5.mat', 'Eb_N0', 'avg_patterns_tested');
semilogy(Eb_N0, avg_patterns_tested, 'bo-', 'linewidth', 1, 'Parent', axes2);
hold on;
%avg_patterns_tested_symbol_GRAND1 = avg_patterns_tested;
clear avg_patterns_tested Eb_N0

xlim(axes2,[20 34]);
ylim(axes2,[10^1 10^5]);
set(axes2,'xtick',20:2:34, 'FontSize',14)
box(axes2,'on');
grid(axes2,'on');

ylabel('Avg. num. of tests per block','FontSize',14, 'Parent', axes2);
xlabel('$E_b/N_0$ (dB)', 'Interpreter','latex','FontSize',16, 'Parent', axes2);

legend1 = legend(axes2, {'RLC & b-GRAND', 'RLC & s-GRAND'});
set(legend1,...
    'Position',[0.195695131057493 0.822832363622651 0.278833587646484 0.132692307692308],...
    'FontSize',12);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.744 0.326923076923077 0.03 0.100000000000015]);

% Create textbox
annotation(figure1,'textbox',...
    [0.622000000000001 0.215153847254245 0.134623291015625 0.10307692197653],...
    'String',{'$w_\mathrm{th}=2$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create arrow
annotation(figure1,'arrow',[0.728000000000002 0.748000000000001],...
    [0.292307692307694 0.342307692307693]);

% Create arrow
annotation(figure1,'arrow',[0.340000000000001 0.374000000000001],...
    [0.415384615384618 0.480769230769233]);

% Create textbox
annotation(figure1,'textbox',...
    [0.242000000000001 0.253615385715785 0.134623291015625 0.10307692197653],...
    'String',{'$w_\mathrm{th}=3$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.240000000000001 0.326692308792708 0.134623291015625 0.10307692197653],...
    'String',{'$w_\mathrm{th}=2$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create arrow
annotation(figure1,'arrow',[0.354 0.400000000000001],...
    [0.323076923076923 0.411538461538464]);

% Create line
annotation(figure1,'line',[0.706000000000001 0.862],...
    [0.745384615384617 0.745384615384617]);

% Create line
annotation(figure1,'line',[0.706000000000001 0.862],...
    [0.676153846153849 0.676153846153849]);

% Create arrow
annotation(figure1,'arrow',[0.85200000000000 0.85200000000000],...
    [0.745384615384617 0.676153846153849],'HeadWidth',5,'HeadStyle','plain',...
    'HeadLength',5);

% Create textbox
annotation(figure1,'textbox',...
    [0.856000000000003 0.630769230769232 0.0819999999999999 0.0720769230769261],...
    'String','3356',...
    'LineStyle','none',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.856000000000003 0.723076923076925 0.082 0.0643846153846176],...
    'String','7654',...
    'LineStyle','none',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.854000000000003 0.680769230769235 0.0819999999999999 0.0643846153846185],...
    'String','-56%',...
    'LineStyle','none',...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

% Create ellipse
annotation(figure1,'ellipse',...
    [0.636000000000001 0.692307692307696 0.03 0.119230769230785]);

% Create arrow
annotation(figure1,'arrow',[0.68 0.652000000000002],...
    [0.861538461538462 0.811538461538467]);

% Create textbox
annotation(figure1,'textbox',...
    [0.680000000000002 0.83438461648502 0.121999999999998 0.088692306591909],...
    'String',{'$w_\mathrm{th}=3$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');
