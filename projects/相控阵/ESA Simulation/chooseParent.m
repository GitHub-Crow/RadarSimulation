function [fa, ma] = chooseParent(currentGeneration, currentGenerationScore, Tournament)
permIndex = randperm(Tournament.TournamentEligible);
competitorFitness = currentGenerationScore(permIndex(1 : Tournament.TournamentCompetitors));
[~, sortedIndex] = sort(competitorFitness);
fa = currentGeneration(:, permIndex(sortedIndex(1)));
ma = currentGeneration(:, permIndex(sortedIndex(2)));
end

