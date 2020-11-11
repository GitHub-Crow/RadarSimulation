delta1 = 0; delta2 = pi*3/4;
SamplePoints = 10000;
AnimatePoints = 5000;
t = linspace (0, 1, SamplePoints);
w = 2*pi*2.8*10^9;
y = cos(w*t + delta1);
z = cos(w*t + delta2);
h1 = animatedline('MaximumNumPoints', AnimatePoints, 'LineStyle', 'none', 'Marker', 'x', 'MarkerFaceColor', 'r');
h2 = animatedline('MaximumNumPoints', 2, 'Marker', 'x', 'MarkerFaceColor', 'r', 'LineWidth', 3, 'Color', 'r');
fig = figure(1);
axis([-1.5 1.5 -1.5 1.5])
axis equal
ST = tic;

for k = 1 : SamplePoints
    addpoints(h1, y(k), z(k));
    addpoints(h2, 0, 0);
    addpoints(h2, y(k), z(k));
    ET = toc(ST);
    if ET > 1/100000
        drawnow 
        ST = tic;
    end
end
