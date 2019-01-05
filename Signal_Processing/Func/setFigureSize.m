function setFigureSize(width, height)

%     h1 = suptitle(fname(9:end));
%     set(h1, 'Interpreter','none','fontsize',4)

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[width height]);
    set(gcf,'PaperPosition',[0 0 width height]);
    box off;