function [ h_fig ] = setFigure(name)
% This function produces an empty figure with noramlized dimensions.

% Parameters
journal_plot_params;

% New figure
h_fig = figure;

height_1_1_bis = height_1_1 - 1;

% Properties
set(h_fig, 'Units', units);
set(h_fig, 'Name', name);
set(h_fig, 'Position', [1, 3, width_1, height_1_1_bis]);
set(h_fig,'Renderer','painters');
set(h_fig,'PaperUnits','centimeters');
set(h_fig,'PaperSize',[width_1 height_1_1_bis]);
set(h_fig,'PaperPositionMode','manual');
set(h_fig,'PaperPosition',[0 0 width_1 height_1_1_bis]);


end

