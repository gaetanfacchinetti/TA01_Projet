clear all;

h_fig = setFigure('Nouvelle Figure')

%Temps  = [32, 15.44, 10.24, 8.04, 6.24, 5.61, 4.87, 5.544, 5.639, 9.229];
%nbProc = [1 3:1:11];
%Rapp   = [0, 71/1360, 162/907, 212/680, 259/544, 289/454, 391/388, 395/340, 433/302, 593/272];

% Rapp   = [71/1360, 162/907, 212/680, 259/544, 289/454, 391/388,...
%     395/340, 433/302, 593/272, 473/243, 521/225, 580/210, 559/194, ...
%     648/181, 630/170, 721/160, 742/151, 766/143];
% nbProc = 3:1:20;
% Temps = [315, 139, 85, 57, 43, 33, 37, 33, 31, 27, 28, 33, 29,...
%     34, 39, 41, 44, 50];

 nbProc = 11:1:20
 Rapp  = [ 301/844, 365/766, 392/703, 355/648, 400/603, 451/563, 417/527, 504/496 548/468, 639/444]
 Temps = [183, 178, 160, 155, 142, 135, 143, 145, 148, 158]

%nbProc = 11:1:20;
%Rapp = [ 82/475, 96/430, 104/392, 97/363, 112/336, 125/313, 118/294, 138/274 146/260];
%Temps = [183, 178, 160 , 152, 142, 135, 143, 145, 148, 158];

[hAx, ht, hr] = plotyy(nbProc, Temps, nbProc, Rapp) ;

X = 0:21;
Y = linspace(max(Temps)/max(Rapp),max(Temps)/max(Rapp),22);

hold on


hp = plot(X,Y,'--');
set(hp, 'Color', 'blue');

set(ht, 'Marker', 'o');
set(ht, 'MarkerEdgeColor', 'red');
set(ht, 'MarkerFaceColor', 'red');
set(hr, 'MarkerSize', 5);
set(ht, 'LineStyle', 'none');

set(hr, 'Marker', 'o');
set(hr, 'MarkerEdgeColor', 'blue');
set(hr, 'MarkerFaceColor', 'blue');
set(hr, 'MarkerSize', 4); 
set(hr, 'LineStyle', 'none');

%set(hAx,'xlim',[0,max(nbProc)+1]);
set(hAx,'xlim',[10,max(nbProc)+1]);
set(hAx(1),'ylim', [130, max(Temps(:))*1.10]);
set(hAx(2),'ylim', [0, max(Rapp(:))*1.10]);
%set(hAx(1), 'YTick', [0,10,20,30]);
set(hAx(1), 'YTick', [0,150, 175, 200]);
%set(hAx(2), 'YTick', [0,1,2]);
set(hAx(2), 'YTick', [0,1,2,3,4,5]);
set(hAx(2), 'YMinorTick', 'off');
set(hAx(2), 'TickLength', [0.001 0.001]);

h_axis = gca;

xlabel(h_axis,'Nombre de processeurs','Interpreter','LaTeX',...
    'FontSize',10);

ylabel(hAx(1),'Temps d''execution (s)','Interpreter','LaTeX',...
    'FontSize',10, 'Color', 'red');
set(hAx(1), 'YColor', 'red');

ylabel(hAx(2),'$\tau$', 'Interpreter', 'LaTeX',...
    'FontSize',10, 'Color', 'blue');
set(hAx(2), 'YColor', 'blue');


journal_plot_params;

set(h_axis, 'Units', 'normalized');
set(h_axis, 'LineWidth', alw);
set(h_axis, 'Position', [left_space+0.03 v_space norm_ax_width-0.01 norm_ax_height]);
set(h_axis, 'FontUnits','points');
set(h_axis, 'FontWeight','normal');
set(h_axis, 'FontSize', fsz);
set(h_axis, 'FontName', f_name);

set(hAx(2), 'Units', 'normalized');
set(hAx(2), 'LineWidth', alw);
set(hAx(2), 'FontUnits','points');
set(hAx(2), 'FontWeight','normal');
set(hAx(2), 'FontSize', fsz);
set(hAx(2), 'FontName', f_name);
