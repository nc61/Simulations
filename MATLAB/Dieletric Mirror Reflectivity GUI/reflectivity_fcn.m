function [wv, R, full_layer_matrix] = reflectivity_fcn(wv_range, npoints, usefile, data)



% Initialize independent variables
wv_start = wv_range(1);
wv_end = wv_range(2);
wv = linspace(wv_start, wv_end, npoints);
k0 = 2*pi./wv;

if (usefile)
    % Open the file and read the data
    fid = fopen(data);
    layers = fscanf(fid, '%f');
    layers = reshape(layers, 9, length(layers)/9)';
    full_layer_matrix = layers;
    fclose(fid);
else
    full_layer_matrix = data;
    layers = data;
end

% Get the substrate and ambient indices of refraction
x_substrate = layers(1, 4);
n_substrate = index2(wv, x_substrate);
n_ambient = layers(end, 6);

% Remove the first (substrate) and last (ambient) layers from the matrix
layers(1,:) = [];
layers(end,:) = [];

% Find the layers and layer groups from the layer matrix
num_unique_layers = size(layers,1);
group_starts = [1; find(1 - diff(layers(:,9)) > 0) + 1];
num_layers_in_groups = [diff(group_starts); num_unique_layers - group_starts(end) + 1];

% Initialize the shape of T
zero_vec = zeros(1, length(k0));
zero_cell = {zero_vec, zero_vec;zero_vec, zero_vec};
T = {zero_cell,zero_cell,zero_cell,zero_cell,zero_cell, zero_cell};

% Initialize T_total
T_total = {1,0;0,1};

% Outer loop through each group
for ind = 1:length(group_starts)
    
    % Find the number of times the layer should repeat
    repetitions = layers(group_starts(ind), 8);
    
    if (ind == 1)
        
        % Create an interface with the substrate
        x1 = layers(1, 4);
        n1 = index2(wv, x1);
        T{ind} = T_interface(n_substrate, n1);
        T{ind} = {1,0;0,1};
    else
        % Create an identity matrix for later multiplication
        T{ind} = {1,0;0,1};
    end
    
    % Inner loop through every layer in group {ind}
    for jnd = 1:num_layers_in_groups(ind)
        
        % Find out which layer (row of matrix) we are in
        layer_index = group_starts(ind) + jnd - 1;
        
        % Compute the index of refraction
        x = layers(layer_index, 4);
        n = index2(wv, x);
        
        % Compute propagation T and multiply it with the overall group T
        len = layers(layer_index, 5)*1e-10;
        T_prop = T_propagation(n, len, k0);
        T{ind} = matmult(T{ind}, T_prop);
        
        % Check if we are in the last layer of group {ind}
        if (jnd == num_layers_in_groups(ind))
            % If we are and the group is repeated, create an interface with
            % the first layer
            if (repetitions > 1)
                first_layer_comp = layers(group_starts(ind), 4);
                T_int = T_interface(n, index2(wv, first_layer_comp));
                T{ind} = matmult(T{ind}, T_int);
                T{ind} = matexp(T{ind}, repetitions);
            end
            
            % If we have another layer after, create the interface.
            % Otherwise, do nothing since we have no defined interface.
            if (group_starts(ind) + jnd < num_unique_layers)
                T_int = T_interface(n, index2(wv, layers(ind + jnd + 1,4)));
                T{ind} = matmult(T{ind}, T_int);
            end
            
        else
            % Otherwise, create an interface with the next layer in the
            % group
            x_next = layers(layer_index + 1, 4);
            n_next = index2(wv, x_next);
            T_int = T_interface(n, n_next);
            T{ind} = matmult(T{ind}, T_int);
            
        end
    end
    
    % Multiply the total transition matrix by the group transition matrix
    T_total = matmult(T_total, T{ind});
end

% Create an interface with the ambient material
x_last = layers(end, 4);
n_last = index2(wv, x_last);
T_int = T_interface(n_last, n_ambient);
T_total = matmult(T_total, T_int);


% Get the reflectivity from the matrix
r = T_total{1,2}./T_total{1,1};
R = abs(r).^2;

end

function T = T_propagation(n, L, k0)
T = {exp(1j*(n.*k0).*L), 0; 0, exp(-1j*(n.*k0).*L)};
end

function T = T_interface(n1, n2)
r = (n2 - n1)./(n1 + n2);
t = sqrt(1 - r.^2);
T = {1./t, r./t; r./t, 1./t};
end

function T = matmult(T1, T2)

T = {0,0;0,0};

T{1,1} = T1{1,1}.*T2{1,1} + T1{1,2}.*T2{2,1};
T{1,2} = T1{1,1}.*T2{1,2} + T1{1,2}.*T2{2,2};
T{2,1} = T1{2,1}.*T2{1,1} + T1{2,2}.*T2{2,1};
T{2,2} = T1{2,1}.*T2{1,2} + T1{2,2}.*T2{2,2};

end

function T_N = matexp(T, N)
T_N = T;

for ind = 1:N-1
    T_N = matmult(T, T_N);
end

end