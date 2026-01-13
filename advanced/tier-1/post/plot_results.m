function plot_results(nodes, elements, U, stress)

scale = 100;
deformed = nodes + scale*[U(1:2:end), U(2:2:end)];

figure;
patch('Faces',elements,'Vertices',deformed,...
      'FaceVertexCData',stress(:,1),...
      'FaceColor','flat');
colorbar;
axis equal;
title('Deformed Shape (σ_x)');
end
