% Run mumps out-of-core (OOC) to reuse factorizations within parfor.
%
% See mumps userguide chpt. 5.20
%
% Issue: MATLAB workers in a parfor loop do not share memory, hence objects
% like decomposition/MUMPS cannot be directly broadcasted to workers
% -> these objects in parts only exist in memory
% Fix: Mumps factorization needs to be stored to disc (serialized) and read
% by the worker (deserialized) using the mumps-own OOC

%% Create test case.

n_sys = 3e3;    % system size
n_A = 3;        % number of different systems
n_RHS = 20;     % number of rhs vectors
n_runs = 5;     % number of run repetitions for timing
n_rep = n_A*30; % number of solving repetitions, each with block-RHS
n_workers = 2;  % number of MATLAB parallel worker

% Create a sparse matrix A with good condition number.
A = cell(n_A, 1);
rng(0815);      % for reproducibility
for ii = 1:length(A)
    A{ii} = sprandn(n_sys, n_sys, 0.2) + speye(n_sys);
    A{ii} = A{ii} * A{ii}';                     % A symmetric
    A{ii} = A{ii} + 10 * speye(size(A{ii}));    % A with improved condition number
end

% Create the RHS matrix B.
B = zeros(n_sys, n_RHS);
for ii = 1:size(B, 2)
    B(randperm(n_sys, 10), ii) = 1;  % B with only 10 entries per column, sparse
end
B = sparse(B);

% Create mapping between Systemmatrices and repetitions.
% -> we want to test, if multiple solver objects could be loaded and used
%    simultaneously within parfor
vwd2rep = reshape(repmat(1:n_A, n_rep / n_A, 1), [], 1);

%% Multiple matrices serial

% Create solver
solver = cell(n_A, 1);
tic;
% Perform decomposition of A once and store on disc
for ii = 1:length(A)
    solver{ii} = MUMPS_OOC(A{ii}, ii);
end
time_deco_store = toc;

% Loop over chunks
X = zeros(n_sys, n_RHS, n_rep);
[time_serial, time_load] = deal(zeros(n_runs, 1));
for j = 1:n_runs
    % Load decompositions from disc
    tic;
    solver = load_all_solver(solver);
    time_load(j) = toc;
    tic;
    for kk = 1:n_rep
        % Solve the system using different precomputed decompositions
        X(:, :, kk) = solver{vwd2rep(kk)}.solve(B);
    end
    time_serial(j) = toc;
end
fprintf(sprintf('2: %.2d s (load + %i FBS) for %i repetitions in serial, multiple matrices\n', ...
        mean(time_load) + mean(time_serial), n_RHS, n_rep));
fprintf(sprintf('deco&store: %.2d   load: %.2d   solve: %.2d\n', ...
        mean(time_deco_store), mean(time_load), mean(time_serial)));

%% Multiple matrices parallel

% Start the parallel pool
par_pool = gcp('nocreate');
if isempty(par_pool)
    parpool('local', n_workers);
elseif par_pool.NumWorkers ~= n_workers
    delete(par_pool);
    parpool('local', n_workers);
end

% % Make sure no relicts from MUMPS itself are transfered
% % Note: Solver object with property ~isempty(solver.id) can not be
%         transfered between workers! Always clean up BEFORE solver
%         variable is defined as parallel.pool.constant!
for ii = 1:length(solver)
    solver{ii}.unload();
end

% Define pool constants that should only be copied once to the worker
worker_B = parallel.pool.Constant(B);
worker_solver = parallel.pool.Constant(@() load_all_solver(solver));

% Loop over chunks
X = zeros(n_sys, n_RHS, n_rep);
[time_parallel] = deal(zeros(n_runs, 1));
for j = 1:n_runs
    tic;
    parfor kk = 1:n_rep
        % Solve the system using different precomputed decompositions
        all_solver = worker_solver.Value;
        local_solver = all_solver{vwd2rep(kk)};
        X(:, :, kk) = local_solver.solve(worker_B.Value);
    end
    time_parallel(j) = toc;
end
fprintf(sprintf('3: %.2d s (load + %i FBS) for %i repetitions in parallel (%d worker), multiple matrices\n', ...
        mean(time_parallel), n_RHS, n_rep, n_workers));

%% Final clean up files

for ii = 1:length(A)
    solver{ii}.delete();
end

%% Helper

function solver = load_all_solver(solver)
    % Restore solver objects on each worker once from !shared! disc

    for ii = 1:length(solver)
        solver{ii}.load();
    end
end
