% Run mumps out-of-core (OOC) to reuse factorizations within parfor.
%
% See mumps userguide chpt. 5.20
% 
% Issue: MATLAB workers in a parfor loop do not share memory, hence objects
% like decomposition/MUMPS cannot be directly broadcasted to workers
% -> these objects in parts only exist in memory
% Fix: Mumps factorization needs to be stored to disc (serialized) and read
% by the worker (deserialized) using the mumps-own OOC

n_workers = 2;
n_rep = 100;   % parallel iterations
n_sys = 3e3;   % system size
n_RHS = 1e1;   % number of rhs vectors
rng(0815);     % for reproducibility

%% Create test case.
fprintf('Define test case and solve\n');

% Create a sparse matrix A with 'good' condition number.
A = sprandn(n_sys, n_sys, 0.2) + speye(n_sys);
A = A * A';                          % A symmetric
A = A + 10 * speye(size(A));         % A with improved condition number

% Create the block-RHS B.
B = zeros(size(A, 2), n_RHS);
for ii = 1:size(B, 2)
    B(randperm(n_sys, 10), ii) = 1;  % B with only 10 entries per column, sparse
end
B = sparse(B);

% Get a reference solution.
X_ref = full(A\B);

%% Create mumps object

% Set directories & prefixes.
save_dir = fullfile(pwd, 'tmp_mumps');  % storage location
save_prefix = 'factors';                % file name
if ~exist(save_dir, 'dir')
    mkdir(save_dir);    
end

% Set required env vars.
% Note: seems to be mandatory, rather than defining them as object
%       parameter.
setenv('MUMPS_SAVE_DIR', save_dir);
setenv('MUMPS_SAVE_PREFIX', save_prefix);

%%  Step 1: factorize and save

% Initialize mumps object.
id = initmumps();      % set default id (JOB = -1), i.e. mumps parameter
id = dmumps(id);       % call wrapper to create mumps object (dmumps = real)
id.SYM = 1;

% Set parameter to enable OOC.
% Note: set BEFORE running mumps with id.JOB = 4! (otherwise high memory required -> why?)
%       set AFTER calling z/dmumps!               (otherwise ICNTL is replaced)
id.ICNTL(22) = 1;  % OOC factorization and solve phases. The complete matrix of factors is written to disk
id.ICNTL(35) = 3;  % In an OOC context, (ICNTL(22)=1) this option enables the user to
                   % write all factors to disk which is not the case with ICNTL(35)=2 since factors in low-rank
                   % form are not written to disk.                                -> required?

% % Set save/restore location as parameter.
% id.SAVE_DIR    = save_dir;
% id.SAVE_PREFIX = save_prefix;

% Factorize OOC.
fprintf('MUMPS: factorize ... ');
tic;
id.JOB = 4;                       % analysis + factorization
id = dmumps(id, A);               % IMPORTANT: pass A for factorization
fprintf(sprintf('done (%.2d s)\n', toc));

% Save the instance to disk.
fprintf('MUMPS: save ... ');
tic;
id.JOB = 7;
id = dmumps(id, A);
fprintf(sprintf('done (%.2d s)\n', toc));

% Keep parameter that are not stored
save_SYM = id.SYM;
save_SYM_PERM = id.SYM_PERM;

% Delete the instance, i.e. free internal memory (saved files remain).
id.JOB = -2;
id = dmumps(id, sparse([], [], [], n_sys, n_sys)); %#ok<*NASGU>

% Sanity check.
fprintf('Saved files in %s:\n', save_dir);
saved = dir(fullfile(save_dir, [save_prefix '*']));
disp({saved.name}');

%%  Step 2: restore and solve (local)

% Initialize.
id = initmumps();
id = dmumps(id);

% % Reset location to above defined paths and solver parameter.
% id.SAVE_DIR    = save_dir;
% id.SAVE_PREFIX = save_prefix;
id.SYM = save_SYM;
id.SYM_PERM = save_SYM_PERM;

% Restore saved instance.
fprintf('MUMPS_ser: Restore ... ');
tic;
id.JOB = 8;
% Note: pass empty sparse of correct size such that the dumumps wrapper
%       don't throw an error (passing A as second argument is mandatory).
id = dmumps(id, sparse([], [], [], n_sys, n_sys));
fprintf(sprintf('done (%.2d s)\n', toc));

% Loop over serial FBS
fprintf(sprintf('MUMPS_ser: solve %i times ... ', n_rep));
tic;
X = zeros(n_sys, n_RHS, n_rep);
for kk = 1:n_rep
    % Set RHS and solve.
    id.RHS = B;
    id.JOB = 3;
    id = dmumps(id, sparse([], [], [], n_sys, n_sys));
    X(:, :, kk) = id.SOL;
end
fprintf(sprintf('done (%.2d s)\n', toc));

% Sanity check.
assert(norm(X(:,:,end) - X_ref) < 1e-13);

% Cleanup object.
id.JOB = -2;
id = dmumps(id, sparse([], [], [], n_sys, n_sys));

%%  Step 2: restore and solve (parallel)

% Start the parallel pool
par_pool = gcp('nocreate');
if isempty(par_pool)
    parpool('local', n_workers);
elseif par_pool.NumWorkers ~= n_workers
    delete(par_pool);
    parpool('local', n_workers);
end

% Set up (large), read-only variables that should not be passed between
% worker local iterations.
% Note: MATLAB's "broadcast warning" means data is copied from the client
%       to every worker FOR EVERY LOOP ITERATION!
worker_B = parallel.pool.Constant(B);

% Loop over parallel FBS
fprintf(sprintf('MUMPS_par: Restore and solve %i times ... ', n_rep));
% Define storage variabe outside parfor, X = id_.SOL; otherwise would only
% exist on workers locally!
X = zeros(n_sys, n_RHS, n_rep);
tic;
parfor kk = 1:n_rep
    % Initialize.
    id_ = initmumps();
    id_ = dmumps(id_);
    
    % % Reset location to above defined paths and solver parameter.
    % id.SAVE_DIR    = save_dir;
    % id.SAVE_PREFIX = save_prefix;
    id_.SYM = save_SYM;
    id_.SYM_PERM = save_SYM_PERM;

    % Restore saved instance.
    id_.JOB = 8;
    % Note: pass empty sparse of correct size such that the dumumps wrapper
    %       don't throw an error (passing A as second argument is mandatory).
    id_ = dmumps(id_, sparse([], [], [], n_sys, n_sys));
    
    % Set RHS and solve.
    id_.RHS = worker_B.Value;
    id_.JOB = 3;
    id_ = dmumps(id_, sparse([], [], [], n_sys,n_sys));
    X(:, :, kk) = id_.SOL;
end
fprintf(sprintf('done (%.2d s)\n', toc));

% Sanity check.
assert(norm(X(:,:,end) - X_ref) < 1e-13);

%% Clean up

% Initialize.
id_ = initmumps();
id_ = dmumps(id_);
id_.SYM = save_SYM;

% Restore saved instance.
% Note: loading before cleanup is mandatory.
id_.JOB = 8;
id_ = dmumps(id_, sparse([], [], [], n_sys, n_sys));

% Remove saved files.
id_.JOB = -3;
id_ = dmumps(id_, sparse([], [], [], n_sys, n_sys));
rmdir(save_dir);

% Cleanup object.
id_.JOB = -2;
id_ = dmumps(id_, sparse([], [], [], n_sys, n_sys));
