function plotMatrix(matrix, txt, opts, coord)
if opts.setFloor == 1
    matrix(matrix < opts.floorValue) = opts.floorValue;
end
figure;clf
set(gcf,'DefaultLineLineWidth',1.5)
if opts.dimension == 2
    imagesc(matrix);hold on
else
    if opts.setAxis == 0
        surf(matrix);hold on
    else
        surf(coord.x, coord.y, matrix);hold on
end
shading interp
set(gca,'FontSize',10,'FontWeight','bold')
title(txt.title);
xlabel(coord.xLabel),ylabel(coord.yLabel);
axis tight
colorbar
if opts.setAxis == 1
    set(gca,'XTick', [coord.xMin: coord.xDelta: coord.xMax])
    set(gca,'YTick', [coord.yMin: coord.yDelta: coord.yMax])
end
set(gcf, 'color', 'white');
end

