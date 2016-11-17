clear all;

h_fig = setFigure('Nouvelle Figure')

Temps  = [15.44, 10.24, 8.04, 6.24, 5.61, 4.87, 5.544, 5.639, 9.229];
nbProc = 3:1:11;
Rapp   = [71/1360, 162/907, 212/680, 259/544, 289/454, 391/388, 395/340, 433/302, 593/272];


[hAx, ht, hr] = plotyy(nbProc, Temps, nbProc, Rapp) ;

X = 0:12;
Y = linspace(max(Temps)/max(Rapp),max(Temps)/max(Rapp),13);

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

%set(hAx,'xlim',[0,max(N(:))+10]);
set(hAx(1),'ylim', [0, max(Temps(:))*1.10]);
set(hAx(2),'ylim', [0, max(Rapp(:))*1.10]);
set(hAx(1), 'YTick', [0,5,10,15]);
set(hAx(2), 'YTick', [0,1,2]);


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
