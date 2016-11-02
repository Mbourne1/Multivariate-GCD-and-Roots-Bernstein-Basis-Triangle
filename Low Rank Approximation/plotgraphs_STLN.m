global SETTINGS

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        fig_name = sprintf('%s : condition',mfilename);
        figure('name',fig_name)
        hold on
        plot(log10(condition))
        xlabel('iteration');
        ylabel('log_{10} condition');
        hold off
        
    case 'n'
        
end

