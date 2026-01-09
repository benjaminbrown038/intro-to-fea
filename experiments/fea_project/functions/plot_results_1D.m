function plot_results_1D(nodes, U, elements, stress)
    figure;
    
    subplot(2,1,1)
    plot(nodes, U, '-o','LineWidth',1.5)
    xlabel('x [m]'); ylabel('Displacement [m]')
    title('Nodal Displacements')
    grid on

    subplot(2,1,2)
    midpoints = (nodes(elements(:,1)) + nodes(elements(:,2)))/2;
    plot(midpoints, stress, '-s','LineWidth',1.5)
    xlabel('x [m]'); ylabel('Stress [Pa]')
    title('Element Stresses')
    grid on
end
