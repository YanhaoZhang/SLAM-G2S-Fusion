function [ ] = Plot_CONS_OBJ( CON, OBJ, Fid )
%--------------------------------------------------------------------------%
%--------------------------------------------------------------------------%
DIM = 1 : size(CON,1);
%-------------------- Line specificaiton for constraints ------------------%
color1 = 'c';
LineSpecCon.LineStyle = '-';
LineSpecCon.LineWidth = 2.0;
LineSpecCon.Color = color1;
LineSpecCon.Marker = 'o';
LineSpecCon.MarkerSize = 4.0;
LineSpecCon.MarkerEdgeColor = color1;
LineSpecCon.MarkerFaceColor = color1;
%-------------------- Line specificaiton for objective --------------------%
color2 = 'r';
LineSpecObj.LineStyle = '-';
LineSpecObj.LineWidth = 8;
LineSpecObj.Color = color2;
LineSpecObj.Marker = 'o';
LineSpecObj.MarkerSize = 8;
LineSpecObj.MarkerEdgeColor = color2;
LineSpecObj.MarkerFaceColor = color2;
%--------------------------------------------------------------------------%
figure(Fid)
%---------------------------- Set Figure Size -----------------------------%
%[a,b,c,d] (a,b) uuper left corner;  c = width, d =height
set(gcf, 'Position', [100, 400, 2100, 400]);


hold on
H2 = plot(DIM, OBJ, LineSpecObj);
H1 = plot(DIM, CON, LineSpecCon);
% H2 = scatter(DIM, OBJ, 140);
% H1 = scatter(DIM, CON, 50, 'filled');
% H2 = stem(DIM, OBJ, 'c', 'filled');
% H1 = stem(DIM, CON, 'm');
hold off
box on
%----------------------------- set axis limits------------------------------%
lim = axis;
lim(1) = 1;  % 1
lim(3) = 0;
lim(2) = size(CON,1);
lim(4) = lim(4);
%--------------%
lim(2) = lim(2) + 0.25;
lim(4) = lim(4)*1.05;
%--------------%
axis(lim);
%--------------------------------------------------------------------------%
%--------------------------------------------------------------------------%
%xlabel('The constraint added for each sub-problem', 'FontSize', 20) % x-axis label
% ylabel('constriant metric/objective growth', 'FontSize', 20) % y-axis label
%--------------------------------------------------------------------------%
% legend([H1, H2], {'Constraint metric','Objective growth'}, ... 
%     'FontSize', 20, ...
%     'Location','northwest', ...
%     'TextColor','blue');
%--------------------------------------------------------------------------%
legend([H1, H2], {'\color{cyan} Constraint metric', ...
                  '\color{red} Objective growth'}, ... 
    'FontSize', 28, ...
    'Location','northwest' ...
    );
%--------------------------------------------------------------------------%
%--------------------------------------------------------------------------%
%              Some very useful setting for figure size                    %    
%--------------------------------------------------------------------------%
%---------------------------- Set yTick Digit -----------------------------%
ytickformat('%05.3f')
%--------------------------------------------------------------------------%
%--------------------------------------------------------------------------%
ax = gca;
ax.YAxis.Exponent = 0;         % remove scientific settings.
ax.FontSize = 28;
ax.LineWidth = 1.5;
%--------------------------------------------------------------------------%
%-------------------------- Set Minimal Margin ----------------------------%
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left+0.005, bottom, ax_width-0.01, ax_height];

tightfig

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];

end

