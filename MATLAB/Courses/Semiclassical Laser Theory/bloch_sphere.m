
% Cut out some unnecessary vectors so it doesn't take so long to plot
sz = round(size(U,2)/50);

% un-normalize so the maximum vector has unity length
U2 = Ww.*U*sqrt(pi)*w_in;
V2 = Ww.*V*sqrt(pi)*w_in;
W2 = Ww.*W*sqrt(pi)*w_in;

% downsample
U2 = U2(:,1:sz:end);
V2 = V2(:,1:sz:end);
W2 = W2(:,1:sz:end);

% Define a zero poitn for drawing the vectors
zero_vec = zeros(1, size(U2,2));

figure(6)

% Play a movie of the bloch sphere as it evolves in time
movie = 0;
if movie == 1
    for ind = 1:length(t_total)
        
        quiver3(zero_vec, zero_vec, zero_vec, U2(ind, :), V2(ind, :), W2(ind, :))
        xlim([-1 1]), xlabel('U');
        ylim([-1 1]), ylabel('V');
        zlim([-1 1]), zlabel('W');
        
        pause(0.01)
    end
else
        % Just plot a single point of interest
        t_plot = 22.06;
    
        ind = max(find(t_total <= t_plot*1e-6));
        quiver3(zero_vec, zero_vec, zero_vec, U2(ind, :), V2(ind, :), W2(ind, :))
        xlim([-1 1]), xlabel('U');
        ylim([-1 1]), ylabel('V');
        zlim([-1 1]), zlabel('W');
end

