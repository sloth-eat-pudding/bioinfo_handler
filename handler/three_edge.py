import itertools

paths = ['rrr', 'rra', 'rar', 'raa', 'arr', 'ara', 'aar', 'aaa']

# Calculate all combinations of 3 distinct paths
path_combinations = list(itertools.combinations(paths, 3))
for i, combo in enumerate(path_combinations):
    print(f"{i+1}: {combo}")