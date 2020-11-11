%% 2D pattern
%% ESA parameters
Array.f = 3;
Array.f0 = 3;
Array.nxele = 2;
Array.nyele = 2;
Array.lambda = 0.3/Array.f;
Array.dx = 1/(1 + sind(50))*Array.lambda;
Array.dy = 1/(1 + sind(50))*Array.lambda;
Array.EF = 1.5;
Array.flag.wgt = 0; % 0 for uniform, 1 for Taylor Weighting
Array.flag.errRand = 1; % 0 for ideal, 1 for random amplitude and phase error
Array.flag.errQuant = 1; % 0 for ideal, 1 for quant
Array.tilt.opt = 2; % 0 for ideal, 1 for roll, 2 for pitch, 3 for yaw
Array.tilt.angle = 60;
Array.Error.amSigma = 0.5;
Array.Error.phSigma = 6;
Array.qBits = 4;

% Taylor weights paramenters
Array.Taylor.nbar = 5;
Array.Taylor.SLL = 30;

% theta angle paramenters
theta.N = 361; 
theta.minAngle = 0; 
theta.maxAngle = 90;
theta.scanAngle = 60;

phi.N = 361;
phi.minAngle = 0;
phi.maxAngle = 360;
phi.scanAngle = 90;

%% computations

if Array.flag.wgt == 0
    Array.amWeights = ones(Array.nyele, Array.nxele);
else
    amWeightsX = getTaylorWeights(Array.nxele, Array.Taylor.SLL, Array.Taylor.nbar);
    amWeightsY = getTaylorWeights(Array.nyele, Array.Taylor.SLL, Array.Taylor.nbar);
    Array.amWeights = amWeightsX*amWeightsY';
end

% set random error or ideal error
if Array.flag.errRand == 1
    Array.Error.amMatdB = Array.Error.amSigma*randn(Array.nyele, Array.nxele);
    Array.Error.amMat = 10.^(Array.Error.amMatdB/10);
    Array.Error.phMat = Array.Error.phSigma*randn(Array.nyele, Array.nxele);
else
    Array.Error.amMat = ones(Array.nyele, Array.nxele);
    Array.Error.amMatdB = 10*log10(Array.Error.amMat);
    Array.Error.phMat = zeros(Array.nyele, Array.nxele);
end
    


theta.vec = linspace(theta.minAngle, theta.maxAngle, theta.N);
phi.vec = linspace(phi.minAngle, phi.maxAngle, phi.N);
[theta.mat, phi.mat] = meshgrid(theta.vec, phi.vec);

% set rotating matrix

if Array.tilt.opt == 0
    Array.rotMat = eye(3);
elseif Array.tilt.opt == 1
    Array.rotMat = [cosd(Array.tilt.angle) -sind(Array.tilt.angle) 0;
                    sind(Array.tilt.angle) cosd(Array.tilt.angle)  0;
                    0                      0                       1];
elseif Array.tilt.opt == 2
    Array.rotMat = [1 0                       0                      ;    
                    0 cosd(Array.tilt.angle)  -sind(Array.tilt.angle);
                    0 -sind(Array.tilt.angle) cosd(Array.tilt.angle) ];
elseif Array.tilt.opt == 3
    Array.rotMat = [cosd(Array.tilt.angle) 0 -sind(Array.tilt.angle);
                    0                      1 0                      ;
                    sind(Array.tilt.angle) 0 cosd(Array.tilt.angle) ;];
end

sine.uMat = sind(theta.mat).*cosd(phi.mat);
sine.vMat = sind(theta.mat).*sind(phi.mat);
sine.wMat = cosd(theta.mat);

sine.u0 = sind(theta.scanAngle)*cosd(phi.scanAngle);
sine.v0 = sind(theta.scanAngle)*sind(phi.scanAngle);
sine.w0 = cosd(theta.scanAngle);

uvwMat = vertcat(sine.uMat(:)', sine.vMat(:)', sine.wMat(:)');
uvwRotMat = Array.rotMat*uvwMat;
sine.uRotMat = reshape(uvwRotMat(1, :), size(sine.uMat));
sine.vRotMat = reshape(uvwRotMat(2, :), size(sine.vMat));
sine.wRotMat = reshape(uvwRotMat(3, :), size(sine.wMat));

% compute AF

if Array.flag.errQuant == 1
    Array.AF = compute2D_AFQuant(Array.amWeights, Array.dx, Array.dy,...
               Array.f, Array.f0, Array.Error.amMat, Array.Error.phMat,...
               sine.uMat, sine.u0, sine.vMat, sine.v0, Array.qBits);
else
    Array.AF = compute2D_AF(Array.amWeights, Array.dx, Array.dy, ...
               Array.f, Array.f0, Array.Error.amMat, Array.Error.phMat,...
               sine.uMat, sine.u0, sine.vMat, sine.v0);
end

[AF.mag, AF.magdB, AF.magdB_Norm] = processMat(Array.AF);
% compute EP
Array.EP = compute2D_EP(theta.mat, Array.EF);
[EP.mag, EP.magdB, EP.magdB_Norm] = processMat(Array.EP);

% compute PAT
Array.PAT = compute2D_PAT(Array.EP, Array.AF);
[PAT.mag, PAT.magdB, PAT.magdB_Norm] = processMat(Array.PAT);

% compute gain
[Array.intGain, Array.idealGain] = compute2D_IntGain(Array.PAT, Array.f, Array.nxele, Array.nyele,...
                                   Array.dx, Array.dy, theta.vec, phi.vec);
[Gain.mag, Gain.magdB, Gain.magdB_Norm] = processMat(Array.intGain);

%% plot
plotCommand.coordType = 0; % 1 for sinespace, 0 for radar coordination

plotCommand.plotAll = 0;
plotCommand.plotError = 0;
plotCommand.plotEP = 0;
plotCommand.plotPAT = 1;
plotCommand.plotAF = 0;
plotCommand.plotGain = 0;

if plotCommand.plotAll || plotCommand.coordType == 0
    coord.xLabel = '\theta_{AZ}';
    coord.yLabel = '\theta_{EL}';
    coord.x = 180/pi*atan2(sind(theta.mat).*cosd(phi.mat), cosd(theta.mat));
    coord.y = 180/pi*asin(sind(theta.mat).*sind(phi.mat));
    coord.xMin = -90;
    coord.xMax = 90;
    coord.yMin = -90;
    coord.yMax = 90;
    coord.xDelta = 30;
    coord.yDelta = 30;
else
    coord.xLabel = 'u';
    coord.yLabel = 'v';
    coord.x = sine.uMat;
    coord.y = sine.vMat;
    coord.xMin = -1;
    coord.xMax = 1;
    coord.yMin = -1;
    coord.yMax = 1;
    coord.xDelta = 0.5;
    coord.yDelta = 0.5;
end
    

if plotCommand.plotAll || plotCommand.plotError == 1
    txt.title = 'Random Amplitude Error(dB)';
    opts.dimension = 2;
    opts.setFloor = 0;
    opts.setAxis = 0;
    plotMatrix(Array.Error.amMatdB, txt, opts, coord);
    
    figure;clf
    binvector=[-2:.1:2];
    set(gcf,'DefaultLineLineWidth',1.5)
    histogram(Array.Error.amMatdB(:),binvector);
    set(gca,'FontSize',14,'FontWeight','bold')
    title('Amplitude Error Distribution')
    xlabel('Amplitude Error Bins (dB)')
    axis tight
    set(gca,'XTick',[-2:0.5:2])
    set(gcf, 'color', 'white');
    
    txt.title = 'Random Phase Error(degree)';
    plotMatrix(Array.Error.phMat, txt, opts)
    
    figure;clf
    binvector=[-20:1:20];
    set(gcf,'DefaultLineLineWidth',1.5)
    histogram(Array.Error.phMat(:),binvector);
    set(gca,'FontSize',14,'FontWeight','bold')
    title('Amplitude Error Distribution')
    xlabel('Amplitude Error Bins (dB)')
    axis tight
    set(gca,'XTick',[-20:1:20])
    set(gcf, 'color', 'white');
end

if plotCommand.plotAll || plotCommand.plotEP == 1
    txt.title = 'Element Pattern';
    opts.dimension = 3;
    opts.setFloor = 1; opts.floorValue = -20;
    opts.setAxis = 1;
    plotMatrix(EP.magdB_Norm, txt, opts, coord);
end

if plotCommand.plotAll || plotCommand.plotAF == 1
    txt.title = 'Array Factor';
    opts.dimension = 3;
    opts.setFloor = 1; opts.floorValue = -50;
    opts.setAxis = 1;
    plotMatrix(AF.magdB_Norm, txt, opts, coord);
end

if plotCommand.plotAll || plotCommand.plotPAT == 1
    txt.title = 'Pattern';
    opts.dimension = 3;
    opts.setFloor = 1; opts.floorValue = -50;
    opts.setAxis = 1;
    plotMatrix(AF.magdB_Norm, txt, opts, coord);
end

if plotCommand.plotAll || plotCommand.plotGain == 1
    txt.title = 'Integrated Gain';
    opts.dimension = 3;
    opts.setFloor = 0; opts.floorValue = -50;
    opts.setAxis = 1;
    plotMatrix(Gain.magdB_Norm, txt, opts, coord);
end
    
