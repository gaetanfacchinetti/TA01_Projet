function journal_axis( h_axis, label_x, label_y )
% This function changes the axes properties to make them suitable for
% article publishing.

journal_plot_params;

l_fsz = 10;

set(h_axis, 'Units', 'normalized');
set(h_axis, 'LineWidth', alw);
set(h_axis, 'Position', [left_space+0.03 v_space norm_ax_width norm_ax_height]);
set(h_axis, 'FontUnits','points');
set(h_axis, 'FontWeight','normal');
set(h_axis, 'FontSize', fsz);
set(h_axis, 'FontName', f_name);


xlabel(h_axis,label_x,'Interpreter','LaTeX',...
    'FontSize',l_fsz);

ylabel(h_axis,label_y,'Interpreter','LaTeX',...
    'FontSize',l_fsz);



end

