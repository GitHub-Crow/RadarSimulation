clc; close all;
%% init parameters

Nele = 100;
EF = 1.2;

% constraints of excitation

magMin = zeros(Nele, 1);
magMax = ones(Nele, 1);
phsMin = zeros(Nele, 1);
phsMax = ones(Nele, 1)*pi;

enforceSymmetry = 1;

% pattern goals

% define upperbound and lowerbound simplely
upperSidelobeGoalU = [-1 -0.03 0.03 0.4 0.4 0.5 0.5 1];
upperSidelobeGoaldB = [-55 -30 -30 -39.5 -60 -60 -42.1 -55];
lowerSidelobeGoalU = [-1 1];
lowerSidelobeGoaldB = [-80 -80];
autoExemptMainBeam = 1; % 1 for exempt main beam from cost function
NPatternPoints = 512;

% optimization parameters

MaxIterations = 5e2;
popSize = 50;
NMarriages = 25;
Tourment.TournamentEligible = 10;
Tourment.TournamentCompetitors = 5;

NOrder = 2*Nele;
Mutation.MutationProb = 0.04;
Mutation.MutationRangeMax = 0.2;
Mutation.MutationRangeDecay = 1.5;
Mutation.NMutation = round(Mutation.MutationProb*NOrder);


NCrossovers = 2;

magRange = magMax - magMin;
phsRange = phsMax - phsMin;
sin_theta = [-1 : 2/NPatternPoints : 1 - 2/NPatternPoints];
cos_theta = sqrt(1 - sin_theta.^2);
% interpolation to the size of theta
lowerbound_dB = interp1(lowerSidelobeGoalU + [0 : size(lowerSidelobeGoalU, 2) - 1]*eps,...
                        lowerSidelobeGoaldB, sin_theta);
upperbound_dB = interp1(upperSidelobeGoalU + [0 : size(upperSidelobeGoalU, 2) - 1]*eps,...
                        upperSidelobeGoaldB, sin_theta);
lowerbound = 10.^(lowerbound_dB/20)';
upperbound = 10.^(upperbound_dB/20)';
eleWgtPower = cos_theta.^EF;
eleWgtVoltage = sqrt(eleWgtPower)';

%% init population

currentGeneration = rand(NOrder, popSize);
currentGenerationScore = zeros(1, popSize);
for kk = 1 : popSize
    magWgts = currentGeneration(1 : Nele, kk).* magRange + magMin;
    phsWgts = currentGeneration(Nele + 1 : end, kk).*phsRange + phsMin;
    wgts = magWgts.*exp(1j*phsWgts);
    currentGenerationScore(kk) = costFunction(wgts, eleWgtVoltage, lowerbound, upperbound,...
                                              autoExemptMainBeam, enforceSymmetry);
end

%% optimization

GenerationAverageSocre = zeros(MaxIterations, 1);
GenerationBestScore = zeros(MaxIterations, 1);
[~, sort_index] = sort(currentGenerationScore);
currentGeneration = currentGeneration(:, sort_index);
currentGenerationScore = currentGenerationScore(:, sort_index);
globalBestScore = currentGenerationScore(1);

for ii = 1 : MaxIterations
    t = (ii - 1)/(MaxIterations - 1);
    Mutation.mutationRange = Mutation.MutationRangeMax * ...
                             (exp(-Mutation.MutationRangeDecay*t) - ...
                             exp(-Mutation.MutationRangeDecay));
    
    % for each marriage
    Children  = zeros(NOrder, NMarriages*2);
    for jj = 1 : NMarriages
        
        % choose parents
        [Dad, Mom] = chooseParent(currentGeneration, currentGenerationScore, Tourment);
        
        % mating
        Son = Dad;
        Daughter = Mom;
        permIndex = randperm(NOrder - 1) + 1;
        crossoverPoints = sort(permIndex(1 : NCrossovers));
        for crossoverCount = 1 : NCrossovers
            transferRange = [crossoverPoints(crossoverCount) : NOrder];
            geneSlice = Son(transferRange);
            Son(transferRange) = Daughter(transferRange);
            Daughter(transferRange) = geneSlice;
        end
        
        % mutate son and daughter
        Children(:, 2 * jj - 1) = mutate(Son, NOrder, Mutation);
        Children(:, 2 * jj) = mutate(Daughter, NOrder, Mutation);
    end
    
    % score children
    ChildrenScores = zeros(1, 2*NMarriages);
    for kk = 1 : 2*NMarriages
        magWgts = Children(1 : Nele, kk).*magRange + magMin;
        phsWgts = Children(Nele + 1 : 2*Nele, kk).*phsRange + phsMin;
        complexWgts = magWgts.*exp(1j*phsWgts);
        ChildrenScores(kk) = costFunction(complexWgts, eleWgtVoltage, ...
                                          lowerbound, upperbound, ...
                                          autoExemptMainBeam, enforceSymmetry);
    end
    
    % survival of the fittset
    
    combinedGeneration = [currentGeneration Children];
    combinedScores = [currentGenerationScore ChildrenScores];
    
    [~, sortIndex] = sort(combinedScores);
    survivalIndex = sortIndex(1 : popSize);  % keep popsize constant
    currentGeneration = combinedGeneration(:, survivalIndex);
    currentGenerationScore = combinedScores(:, survivalIndex);
    globalBestScore = currentGenerationScore(1);
    
    % save performance
    GenerationAverageSocre(ii) = mean(currentGenerationScore);
    GenerationBestScore(ii) = globalBestScore;
end

GlobalBestIndividual = currentGeneration(:, 1);

%% calculate best beamforming pattern

NPatternPointsFine = NPatternPoints*8;
sine_theta_fine = [-1 : 2/NPatternPointsFine : 1 - 2/NPatternPointsFine];
cosine_theta_fine = sqrt(1 - sine_theta_fine.^2);
eleWgtPowerFine = cosine_theta_fine.^EF;
eleWgtsVoltageFine = sqrt(eleWgtPowerFine);

complexWgts = (GlobalBestIndividual(1 : Nele).*magRange + magMin).*...
              exp(1j*(GlobalBestIndividual(Nele + 1 : 2*Nele).*phsRange + phsMin));

if enforceSymmetry
    complexWgts(Nele/2 + 1 : Nele) = flipud(complexWgts(1 : Nele/2));
end

complexWgts = complexWgts/max(abs(complexWgts));
BestPattern = eleWgtsVoltageFine.*fftshift(fft(complexWgts, NPatternPointsFine));
BestPattern_dB = 20*log10(abs(BestPattern) + eps);
BestPatternNorm_dB = BestPattern_dB - max(BestPattern_dB);

%% plot global best pattern

figure ;
fontsize = 12;
subplot(2, 2, 1)
set(gca, 'fontsize', fontsize);
plot([1 : MaxIterations], GenerationBestScore,...
     [1 : MaxIterations], GenerationAverageSocre,...
     'linewidth', 2);
xlabel('Iteraion #');
ylabel('cost function');
legend('Best', 'Average');
axis tight;

subplot(2, 2, 2)
eles = linspace(1, Nele, Nele);
yyaxis left
plot(eles, abs(complexWgts), 'Marker', '.');
ylabel('Amplitude', 'fontsize', fontsize);
yyaxis right
plot(eles, angle(complexWgts) - min(angle(complexWgts)), 'Marker', '+');
ylabel('Phase', 'fontsize', fontsize);
title('Aperture Weights from GA', 'fontsize', fontsize);

subplot(2, 1, 2)
set(gca, 'fontsize', fontsize);
plot(sine_theta_fine, BestPatternNorm_dB, 'LineWidth', 2); hold on
plot(lowerSidelobeGoalU, lowerSidelobeGoaldB, 'r-.', 'LineWidth', 2);
plot(upperSidelobeGoalU, upperSidelobeGoaldB, 'r:', 'LineWidth', 2);
axis tight
ylim([-80 0]);
title('Far Field Pattern');
ylabel('Magnitude(dB)');
xlabel('sin(\theta)');