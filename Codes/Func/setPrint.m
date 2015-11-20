function setPrint(width, height, fname, printType)

%     h1 = suptitle(fname(9:end));
%     set(h1, 'Interpreter','none','fontsize',4)

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[width height]);
    set(gcf,'PaperPosition',[0 0 width height]);
    box off;
    
    switch upper(printType)
        case 'PNG'
            print('-dpng',[fname '.png'])
        case 'PDF'
            print('-dpdf',[fname '.pdf'])
        case 'TIF'
            print('-dtiff','-r600',[fname '.tif'])
        case 'EPS'
            print('-depsc',[fname '.eps'])
%         case 'FIG'
%             print('-dtiff','-r600',[fname '.tif'])
    end
    
    
    