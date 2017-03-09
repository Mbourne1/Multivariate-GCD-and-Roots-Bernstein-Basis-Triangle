global SETTINGS

if( SETTINGS.PLOT_GRAPHS)
    
    
    fig_name = sprintf('%s : condition',mfilename);
    figure('name',fig_name)
    hold on
    plot(log10(condition),'-s','DisplayName','Condition')
    xlabel('iteration');
    ylabel('log_{10} condition');
    hold off
    
    
    
end

