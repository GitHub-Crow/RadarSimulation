function output = mutate(input, NOrder, Mutation)
permIndex = randperm(NOrder);
mutationDecision = zeros(NOrder, 1);
mutationDecision(permIndex(1 : Mutation.NMutation)) = 1;
mutatedMax = min(input + Mutation.mutationRange/2, 1);
mutatedMin = max(input - Mutation.mutationRange/2, 0);
mutatedGene = rand(NOrder, 1).*(mutatedMax - mutatedMin) + mutatedMin;
output = input.*(1 - mutationDecision) + mutatedGene.*mutationDecision;
end

