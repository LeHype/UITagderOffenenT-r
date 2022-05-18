function weighted_stages_sum = compWeightedStagesSum(stages,s,b)

weighted_stages_sum = 0;
for ii=1:s
    weighted_stages_sum = weighted_stages_sum + b(ii)*stages(:,ii);
end

end