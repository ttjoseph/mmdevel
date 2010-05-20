function [avg, sd] = rescorrel_similarity_expected_value(a, b)

num_trials = 1000;
vals = zeros(num_trials, 1);

for i = 1:num_trials
    vals(i) = rescorrel_similarity(randomize_symmetric_matrix(a), randomize_symmetric_matrix(b));
end

avg = mean(vals);
sd = std(vals);