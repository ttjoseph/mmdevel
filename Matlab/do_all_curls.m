% Runs curl calcuation for several essential modes
subplot_rows = 5;
subplot_cols = 4;
subplot_i = 1;
figure;
for this_ev = 1:10
    the_curl = do_hayward_curl(this_ev, ev_N, ev_C, ev_CA, ev_CB, solute_N, solute_C, solute_CA, solute_CB);
    filename = ['curl.ev', num2str(this_ev), '.out'];
    save(filename, 'the_curl', '-ascii');
    % Create a similarity matrix
    num_atoms = size(the_curl, 1); % Redundant but whatever
    eucdist = zeros(num_atoms, num_atoms); % Euclidean distance (pairwise)
    curlsim = zeros(num_atoms, num_atoms); % Curl similarity (pairwise)
    evsim = zeros(num_atoms, num_atoms); % EV component vector similarity (pairwise)
    the_ev = reshape(ev_CA(:, this_ev), 3, num_atoms);
    the_ev = the_ev';
    for i = 1:num_atoms
        for j = 1:i
            eucdist(i, j) = sqrt(sum((solute_CA(i, :) - solute_CA(j, :)) .^ 2));
            eucdist(j, i) = eucdist(i, j);
            a = the_curl(i, :) / norm(the_curl(i, :));
            b = the_curl(j, :) / norm(the_curl(j, :));
            curlsim(i, j) = dot(a, b);
            curlsim(j, i) = curlsim(i, j);

            a = the_ev(i, :) / norm(the_ev(i, :));
            b = the_ev(j, :) / norm(the_ev(j, :));
            evsim(i, j) = dot(a, b);
            evsim(j, i) = evsim(i, j);
        end
    end
    
    % Calculate the dissimilarity matrix!
    % This is not that important!!
    eucdist = eucdist / max(max(eucdist));
    %similarity = (1 ./ eucdist) .* curlsim;
    % dissimilarity = (eucdist > 0.25) .* (curlsim < 0.5);
    % dissimilarity = eucdist .* curlsim;
    %similarity(isnan(similarity)) = 0;
    %similarity(isinf(similarity)) = 0;
    % dissimilarity(isnan(dissimilarity)) = 0;
    % dissimilarity(isinf(dissimilarity)) = 0;

    filename = ['evsim.ev', num2str(this_ev), '.out'];
    save(filename, 'evsim', '-ascii');
    % filename = ['dissimilarity.ev', num2str(this_ev), '.out'];
    % save(filename, 'dissimilarity', '-ascii');
    filename = ['curlsim.ev', num2str(this_ev), '.out'];
    save(filename, 'curlsim', '-ascii');
    
    % Save the eigenvectors themselves in a format friendly to
    % the VMD arrow drawing script generator
    ev = reshape(ev_CA(:, this_ev), 3, num_atoms);
    ev = ev';
    save(['ev', num2str(this_ev), '.dat'], 'ev', '-ascii');
    
    % Plot stuff so we can save it later
    subplot(subplot_rows, subplot_cols, subplot_i);
    imagesc(evsim);
    ylabel(['EV ', num2str(this_ev)]);
    subplot_i = subplot_i + 1;
    subplot(subplot_rows, subplot_cols, subplot_i);
    imagesc(curlsim);
    subplot_i = subplot_i + 1;
end
me = char(pwd());
me = me(length(me)-2:length(me));
suplabel([me, ' EV and curl similarity matrices'], 't');
fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [8.5 11]);
print('-dpdf', ['sim_', me, '.pdf']);
