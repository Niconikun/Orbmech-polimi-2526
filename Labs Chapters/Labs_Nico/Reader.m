fprintf('\nFirst few time points:\n');
for i = 1:min(5, length(time_utc))
    fprintf('Point %d: %s\n', i, char(time_utc(i)));
end

fprintf('\nKeplerian elements matrix size: %d x %d\n', size(keplerian_elements));
