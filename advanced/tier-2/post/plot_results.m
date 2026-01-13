function plot_results(nodes, elements, U, stress)

scale = 100;
def = nodes + scale*[U(1:2:end), U(2:2:end)];

figure;
patch('Faces',elements,'Vertices',def,...
      'FaceVertexCData',stress(:,1),...
      'FaceColor','flat','EdgeColor','k');
axis equal; colorbar;
title('Deformed Shape with \sigma_x');
end
