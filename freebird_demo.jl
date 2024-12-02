### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ ad6b7485-5af8-4cd3-867a-dc3e53cf9cc0
begin
	# begin
	#     import Pkg
	#     # careful: this is _not_ a reproducible environment
	#     # activate the global environment
	#     Pkg.activate("/Users/myang/Work/git/Demo-FreeBird/Demo")
	# using Plots, PlutoUI, LinearAlgebra
	# using Unitful
	
	# import PlutoUI: combine
	
	# plotly()
	
	# end
	import Pkg; Pkg.add(url="https://github.com/wexlergroup/FreeBird.jl",rev="feature/analysis-and-debug")
	using Plots, PlutoUI, LinearAlgebra
	using Unitful
	using FreeBird

	plotly()
end

# ╔═╡ 82211f2c-c971-48f8-a7a3-e980da7e53c0
html"""<style>
main {
    max-width: 900px;
	font-size: large;
}
"""

# ╔═╡ 4e860050-5b3a-494a-8ce8-e7b5a3ebedf2
md"""
# FreeBird.jl Demo
"""

# ╔═╡ e32b01e1-9c4d-4dd9-85e9-286fd16ed0fa
md"""
# Outline
- Lennard-Jones Potential (the energy calculator)
- What is a walker? (the strcutures/configurations)
- Monte Carlo and Nested sampling
"""

# ╔═╡ 42a7c952-e8e5-43ff-8247-520d6fad4be8
md"""
## Lennard-Jones Potential
"""

# ╔═╡ 8941699d-714b-4072-aadd-f15fb2e4adb4
md"""
$V_\mathrm{LJ}(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]$
"""

# ╔═╡ c61a733f-8a45-4481-a866-7b07d069cfd3
md"""
- ``\epsilon`` = $(@bind epsilon Slider(1:0.1:15)), defines how deep your potential well is
- ``\sigma`` = $(@bind sigma Slider(1:0.1:4)), at which the potential energy is zero ('size of the particle')
- `cutoff` = $(@bind cutoff Slider(1:0.1:4,default=4)) ``\sigma``, the potential energy is set to zero after the cutoff
- `shift`? $(@bind shift CheckBox()) If true, the whole potential will be shifted at the cutoff
"""

# ╔═╡ 3187ac88-e235-48fb-ad52-12b3fd7709ea
lj_potential = LJParameters(epsilon=epsilon,sigma=sigma,cutoff=cutoff,shift=shift)

# ╔═╡ ff1cf1c9-c2c9-4218-8b9d-a8c826c7d898
plot((0.8:0.01:4)*u"Å", r -> lj_energy(r,lj_potential), ylimits=(-epsilon-0.5,2))

# ╔═╡ e8ad89be-bbc3-4c24-800c-b3dd1400c82d
md"""
```julia
struct LJParameters
    epsilon::typeof(1.0u"eV")
    sigma::typeof(1.0u"Å")
    cutoff::Float64
    shift::typeof(0.0u"eV")
end
```
"""

# ╔═╡ 25148d81-f849-4379-bb47-a0f26af0e257
md"## Let's build some walkers"

# ╔═╡ fe5a1609-c37e-4d75-bbe5-813454091d83
md"""
But before that, let's define a new set of Lennard-Jones parameters.
"""

# ╔═╡ 2216fc0e-d8c3-4cfa-a003-f54701fcea66
lj = LJParameters()

# ╔═╡ 161a24f5-e5f1-4174-8236-ccb57d75063c
md"""
Why does this line work?
"""

# ╔═╡ 257840d3-1004-418b-b732-200bc65fb3a0
md"""
```julia
function LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
    if isa(shift,Bool) && shift
        shiftenergy = lj_energy(Float64(epsilon)*u"eV", Float64(sigma)*u"Å", cutoff*Float64(sigma)*u"Å")  
    else
        shiftenergy = Float64(shift)*u"eV"
    end
    return LJParameters(epsilon*u"eV", sigma*u"Å", cutoff, shiftenergy)
end
```
"""

# ╔═╡ 0d37c8c6-f858-4e58-bf41-aae0a6bee7ad
md"""
## Now let's define some walkers
"""

# ╔═╡ c74d16b3-8424-4ba3-bb0c-d7434cabddcc
md"""
First, let's generate some random initial configurations, with 6 hydrogen atoms in a box:
"""

# ╔═╡ 9e91b9d0-6f3c-4a12-8bab-9c8b6f28857e
hydrogens = generate_initial_configs(5::Int64, 100.0::Float64, 6::Int64; particle_type=:H) 

# ╔═╡ de82ef79-26da-4863-bbee-82e65d01be80
md"""
They are configurations, not walkers yet...
"""

# ╔═╡ 391aa43e-a4b3-4e23-a4eb-55c34af2cd51
md"""
## Now let's define some walkers
"""

# ╔═╡ 96978dc5-fa1e-4947-8cab-b969a67468a6
md"""
Warpping the configurations into the `AtomWalker` objects:
"""

# ╔═╡ 205e7cf1-c435-4433-8f9a-9f833197b5e4
walkers = [AtomWalker(config) for config in hydrogens] # comprehension

# ╔═╡ 552a77eb-f578-4497-a580-674ad9702d44
md"""
## Define a liveset
"""

# ╔═╡ ca3cd8bd-3057-4624-b834-35d4e76545b6
md"""
We have previosuly defined a set of walkers, and a set of Lennard-Jones parameters. Let's combine them to make a `liveset`.
"""

# ╔═╡ 3ff30da5-b855-4b74-9d6f-81c7f5f610c4
liveset = LJAtomWalkers(walkers,lj)

# ╔═╡ 9418bea4-cd5b-47eb-97fd-d3e03ced0a8c
md"""
Now we have everything needed for a nested sampling calculation!
"""

# ╔═╡ 30095fa5-cc4b-4260-892e-f5defa326a3c
md"""
## Nested Sampling!
"""

# ╔═╡ eed6a262-f804-419d-ba96-9db596202508
mc = MCRandomWalkClone()

# ╔═╡ 8a2920e2-9058-4e29-8560-8a249530178f
ns_params=NestedSamplingParameters(200::Int64, 0.1::Float64, 0.01::Float64, 1e-5::Float64, 1.0::Float64, 0::Int64, 200::Int64)

# ╔═╡ f8d5538f-0621-4c11-875d-729af53123da
n_iters = 20_000

# ╔═╡ 6999a230-82d1-480d-a0fc-a7bdada6fdb6
save = SaveEveryN(n=20_000) 

# ╔═╡ 354c14c8-df26-460e-bb9b-c3945561d81d
md"""
## Run the code
"""

# ╔═╡ 966a5497-2b72-49bf-9b76-20279a795d22
new_walkers = AtomWalker.(generate_initial_configs(120, 100.0, 6))

# ╔═╡ 68f4f490-f80a-4269-9000-e43ae269544a
begin
	if isfile("output_df.csv")
		run(`rm output_df.csv output.ls.extxyz output.traj.extxyz`) # delete previous output
	end
	new_ls = LJAtomWalkers(new_walkers, lj)
	energies, live, _ = nested_sampling_loop!(new_ls, ns_params, n_iters, mc, save)
end

# ╔═╡ 1beb8de4-bf20-44d3-8ff8-19a5dcfa6f94
# energies

# ╔═╡ f2394846-280f-4c1e-9a82-5c78711452d7
new_ls.walkers

# ╔═╡ fa2a0771-d5b5-44fd-aa50-69a185377aaf
md"""
## `pymatnest` vs. `FreeBird.jl`: different design philosophy
"""

# ╔═╡ c511a108-5b7c-4b45-8789-d6f6d55d583e
md"""
Main nested sampling loop in `pymatnest`: 741 lines
"""

# ╔═╡ 27f32bc0-b10c-4d33-b5d6-46393e3f1a24
md"""
```Python
def do_ns_loop():
    '''
    This is the main nested sampling loop, doing the iterations.
    '''
    global print_prefix
    global cur_config_ind

    # cant keep atoms fixed and change the simulation cell at the moment
    if (movement_args['keep_atoms_fixed'] > 0
            and movement_args['n_cell_volume_steps'] > 0):
        exit_error("cant keep atoms fixed and change the simulation cell\n", 11)

    nD = 3  # dimensionality of the system
    if movement_args['2D']:  # simulate a 2D system only
        nD = 2
    if rank == 0:
        if energy_io.tell() == 0:  # tell returns the current stream position
            if movement_args['do_velocities']:
                nExtraDOF = 0
            else:
                if movement_args['keep_atoms_fixed'] > 0:
                    nExtraDOF = (n_atoms-movement_args['keep_atoms_fixed'])*nD
                else:
                    nExtraDOF = n_atoms*nD
            energy_io.write("%d %d %d %s %d\n" % (ns_args['n_walkers'], ns_args['n_cull'], nExtraDOF, movement_args['MC_cell_flat_V_prior'], n_atoms))

    ## print(print_prefix, ": random state ", np.random.get_state())
    ## if rank == 0:
        ## print(print_prefix, ": common random state ", common_random_state)

    if ns_args['debug'] >= 10 and size <= 1:
        for at in walkers:
            at.info['n_walks'] = 0

    for at in walkers:
        at.info['KEmax'] = KEmax
        if movement_args['MC_cell_P'] > 0:
            print(rank, ": initial enthalpy ", at.info['ns_energy'], " PE ", eval_energy_PE(at), " KE ", eval_energy_KE(at), " PV ", eval_energy_PV(at), " mu ", eval_energy_mu(at), " vol ", at.get_volume())
        else:
            print(rank, ": initial enthalpy ", at.info['ns_energy'], " PE ", eval_energy_PE(at), " KE ", eval_energy_KE(at), " mu ", eval_energy_mu(at), " vol ",at.get_volume())
    sys.stdout.flush()

    # stats for purpose of adjusting step size
    walk_stats_adjust = {}
    # stats for purpose of monitoring acceptance rates
    walk_stats_monitor = {}
    zero_stats(walk_stats_adjust, movement_args)
    zero_stats(walk_stats_monitor, movement_args)

    initial_time = time.time()
    prev_time = initial_time
    step_size_setting_duration = 0.0
    total_step_size_setting_duration = 0.0

    Emax_of_step = None
    Emax_save = []
    i_ns_step_save = []
    traj_walker_list = []
    E_dump_list = []
    E_dump_list_times = []

    verbose = False

    # to avoid errors of unassigned values, if in case of a restart the final number of iter is the same as the starting, stop.
    if start_first_iter == ns_args['n_iter']:
        print("WARNING: Increase the n_iter_times_fraction_killed variable in the input if you want NS cycles to be performed.")
        exit_error("starting iteration and the total number of required iterations are the same,hence no NS cycles will be performed\n", 11)

    last_log_X_n = 0.0
    i_range_mod_n_cull = np.array(range(ns_args['n_cull']))
    i_range_plus_1_mod_n_cull = np.mod(np.array(range(ns_args['n_cull']))+1, ns_args['n_cull'])
    log_X_n_term = np.log(ns_args['n_walkers']-i_range_mod_n_cull) - np.log(ns_args['n_walkers']+1-i_range_mod_n_cull)
    log_X_n_term_cumsum = np.cumsum(log_X_n_term)
    log_X_n_term_cumsum_modified = log_X_n_term_cumsum - np.log(ns_args['n_walkers']+1-i_range_plus_1_mod_n_cull)
    log_X_n_term_sum = log_X_n_term_cumsum[-1]
    if ns_args['converge_down_to_T'] > 0:
        converge_down_to_beta = 1.0/(ns_args['kB']*ns_args['converge_down_to_T'])
        log_Z_term_max = np.NINF

    prev_snapshot_iter = None
    pprev_snapshot_iter = None
    last_snapshot_time = time.time()

    # for estimating current temperature from d log Omega / d E
    if ns_args['T_estimate_finite_diff_lag'] > 0:
        log_alpha = np.log(float(ns_args['n_walkers']+1-ns_args['n_cull'])/float(ns_args['n_walkers']+1))
        Emax_history = collections.deque(maxlen=ns_args['T_estimate_finite_diff_lag'])

    if ns_analyzers is not None:
        for (ns_analyzer, ns_analyzer_interval) in ns_analyzers:
            ns_analyzer.analyze(walkers, -1, "NS_loop start")

    # START MAIN LOOP
    i_ns_step = start_first_iter
    while ns_args['n_iter'] < 0 or i_ns_step < ns_args['n_iter']:

        if ns_args['debug'] == -5:
            print(i_ns_step, rank, " ".join(["{:.2f}".format(eval_energy(x))
                                             for x in walkers]))

        check_memory.check_memory("start_ns_main_loop")
        print_prefix = "%d NS %d" % (rank, i_ns_step)

        if ns_args['debug'] >= 4 and ns_args['track_configs']:
            for at in walkers:
                print(print_prefix, "INFO: 10 config_ind ", at.info['config_ind'], " from ", at.info['from_config_ind'], " at ", at.info['config_ind_time'])

        if movement_args['adjust_step_interval'] < 0:
            zero_stats(walk_stats_adjust, movement_args)
        if movement_args['monitor_step_interval'] < 0:
            zero_stats(walk_stats_monitor, movement_args)

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE START 00 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE START 01 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X START 02 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        # get list of highest energy configs
        (Emax, Vmax, cull_rank, cull_ind) = max_energy(walkers, n_cull)
        Emax_next = Emax[-1]
        if rank == 0 and Emax_of_step is not None and Emax[0] > Emax_of_step:
            print(print_prefix, ": WARNING: energy above Emax ", Emax_of_step, " bad energies: ", Emax[np.where(Emax > Emax_of_step)], cull_rank[np.where(Emax > Emax_of_step)], cull_ind[np.where(Emax > Emax_of_step)])
            # comm.barrier()
            # exit_error("Energy above Emax\n", 5)

        if rank == 0 and (i_ns_step > start_first_iter and Emax_next >= Emax_of_step):
            print("WARNING: Emax not decreasing ", Emax_of_step, Emax_next)
        Emax_of_step = Emax_next

        if ns_args['min_Emax'] is not None and Emax_of_step < ns_args['min_Emax']:
            if rank == 0:
                # if the termination was set by a minimum energy, and it is reached, stop.
                print("Leaving loop because Emax=", Emax_of_step, " < min_Emax =", ns_args['min_Emax'])
            i_ns_step += 1  # add one so outside loop when one is subtracted to get real last iteration it's still correct
            break

        if rank == 0:
            cur_time = time.time()
            output_this_iter = (cur_time > prev_time+60 or i_ns_step == 0 or i_ns_step == ns_args['n_iter'] or (ns_args['n_iter'] > 0 and i_ns_step % max(int(ns_args['n_iter']/1000), 1) == 0))
        else:
            output_this_iter = False

        if ns_args['converge_down_to_T'] > 0:
            # see ns_analyse.py calc_log_a() for math
            log_a = log_X_n_term_sum*i_ns_step + log_X_n_term_cumsum_modified
            # DEBUG if rank == 0:
                # DEBUG for ii in range(len(log_a)):
                    # DEBUG print i_ns_step, "log_a, beta, Es, beta*Es ", log_a[ii], beta, Emax[ii], beta*Emax[ii]
            log_Z_term_max = max(log_Z_term_max, np.amax(log_a - converge_down_to_beta * Emax))
            log_Z_term_last = log_a[-1]-converge_down_to_beta*Emax[-1]
            if output_this_iter:
                print("log_Z_term max ", log_Z_term_max, "last ", log_Z_term_last, "diff ", log_Z_term_max-log_Z_term_last)
            if log_Z_term_last < log_Z_term_max - 10.0:
                if rank == 0:
                    print(print_prefix, "Leaving loop because Z(%f) is converged" % ns_args['converge_down_to_T'])
                i_ns_step += 1  # add one so outside loop when one is subtracted to get real last iteration it's still correct
                break

        if ns_args['T_estimate_finite_diff_lag'] > 0:
            Emax_history.append(Emax_of_step)
        if output_this_iter:
            if ns_args['T_estimate_finite_diff_lag'] > 0 and len(Emax_history) > 1:
                beta_estimate = (len(Emax_history)-1)*log_alpha/(Emax_history[-1]-Emax_history[0])
                T_estimate = 1.0/(ns_args['kB']*beta_estimate)
            else:
                T_estimate = -1
            print(i_ns_step, "Emax_of_step ", Emax_of_step, "T_estimate ", T_estimate, " loop time ", cur_time-prev_time-step_size_setting_duration, " time spent setting step sizes: ", step_size_setting_duration)
            sys.stdout.flush()
            prev_time = cur_time
            step_size_setting_duration = 0.0

        entries_for_this_rank = np.where(cull_rank == rank)[0]
        cull_list = cull_ind[entries_for_this_rank]
        if rank == 0 and ns_args['debug'] >= 4 and len(cull_ind[entries_for_this_rank]) > 0:
            print(print_prefix, "INFO: 20 cull ", cull_ind[entries_for_this_rank], " on ", rank)

        # record Emax walkers energies
        if rank == 0:
            for (E, V) in zip(Emax, Vmax):
                energy_io.write("%d %.60f %.60f\n" % (i_ns_step, E, V))
            energy_io.flush()

            ## Save the energies and corresponding iteration numbers in a list then print(them out only when printing a snapshot)
            #Emax_save.extend(Emax)
            #i_ns_step_save.extend(n_cull*[i_ns_step])
            ## if it is time to print((i.e. at the same iteration when a snapshot is written, or at every iter if no snapshots - for smooth restarts))
            #if ns_args['snapshot_interval'] < 0 or i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1:
            #    for istep,E in zip(i_ns_step_save,Emax_save):
            #        energy_io.write("%d %.60f\n" % (istep, E))
            #    energy_io.flush()
            #    #empty the save lists, so they are ready for the next bunch of saved energies
            #    Emax_save[:]=[]
            #    i_ns_step_save[:]=[]

        # record Emax walkers configurations
        if cull_list is not None:
            for (i, global_n_offset) in zip(cull_list, entries_for_this_rank):
                if ns_args['debug'] >= 10 and size <= 1:
                    print(print_prefix, "walker killed at age ", walkers[i].info['n_walks'])
                # store culled config in list to be written (when
                # snapshot_interval has passed) every traj_interval steps
                global_n = i_ns_step*n_cull + global_n_offset
                if ns_args['traj_interval'] > 0 and global_n % ns_args['traj_interval'] == ns_args['traj_interval']-1:
                    walker_copy = walkers[i].copy()
                    walker_copy.info['volume'] = walker_copy.get_volume()
                    walker_copy.info['ns_P'] = movement_args['MC_cell_P']
                    walker_copy.info['iter'] = i_ns_step
                    walker_copy.info['config_n_global'] = global_n
                    if walker_copy.has('masses') and walker_copy.has('momenta'):
                        walker_copy.info['ns_KE'] = walker_copy.get_kinetic_energy()

                    traj_walker_list.append(walker_copy)

                # if tracking all configs, save this one that has been culled
                if track_traj_io is not None:
                    at = walkers[i].copy()
                    at.info['culled'] = True
                    ase.io.write(track_traj_io, at, parallel=False, format=ns_args['config_file_format'])

        if ns_args['E_dump_interval'] > 0 and i_ns_step % ns_args['E_dump_interval'] == 0:  # ns_args['E_dump_interval']-1:
            if walkers[0].has('masses') and walkers[0].has('momenta'):
                E_dump_list.append([w.info['ns_energy'] - w.get_kinetic_energy() for w in walkers])
            else:
                E_dump_list.append([w.info['ns_energy'] for w in walkers])
            E_dump_list_times.append(i_ns_step)

        if ns_args['traj_interval'] > 0:
            for at in traj_walker_list:
                ase.io.write(traj_io, at, parallel=False, format=ns_args['config_file_format'])
            traj_io.flush()
            traj_walker_list=[]

        # print(the recorded Emax walkers configurations to output file)
        if (ns_args['snapshot_interval'] < 0 or i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1 or
            (ns_args['snapshot_seq_pairs'] and i_ns_step > 0 and i_ns_step%ns_args['snapshot_interval'] == 0) ) :
            ##NB if ns_args['traj_interval'] > 0:
                ##NB for at in traj_walker_list:
                    ##NB ase.io.write(traj_io, at, parallel=False, format=ns_args['config_file_format'])
                ##NB traj_io.flush()
                ##NB traj_walker_list=[]
            if ns_args['E_dump_interval'] > 0:
                if comm is not None:
                    E_dump_list_all = np.array(comm.allgather(E_dump_list))
                else:
                    E_dump_list_all = np.array(E_dump_list)
                if rank == 0:
                    for i in range(E_dump_list_all.shape[1]):
                        E_dump_io.write("step %d\n" % E_dump_list_times[i])
                        if len(E_dump_list_all.shape) == 3:
                            np.savetxt(E_dump_io, E_dump_list_all[:,i,:])
                        else:
                            np.savetxt(E_dump_io, E_dump_list_all[i,:])
                    E_dump_io.flush()
                E_dump_list = []
                E_dump_list_all = None
                E_dump_list_times = []

        # calculate how many will be culled on each rank
        n_cull_of_rank = np.array([sum(cull_rank == r) for r in range(size)])

        # label configs to be culled
        status = np.empty((size, n_walkers), np.object_)
        status[:, :] = ''
        for r in range(size):
            status[r, cull_ind[np.where(cull_rank == r)[0]]] = 'c_t'

        if ns_args['debug'] >= 10:
            initial_PE_loc = [eval_energy(at, do_KE=False) for at in walkers]
            initial_E_loc = [eval_energy(at) for at in walkers]
            if comm is not None:
                initial_PE = np.array(comm.allgather(initial_PE_loc)).flatten()
                initial_E = np.array(comm.allgather(initial_E_loc)).flatten()
            else:
                initial_PE = np.array(initial_PE_loc)
                initial_E = np.array(initial_E_loc)
            initial_changed = initial_PE[np.where(status.flatten() == 'c_t')]
            initial_unchanged = initial_PE[np.where(status.flatten() == '')]

        if ns_args['debug'] >= 30:
            for r in range(len(status)):
                print(print_prefix, ": initial status ", r,
                      [s for s in status[r, :]])

        # find load balance by cloning on top of excess maxima
        recv_ind = []
        recv_rank = []
        send_ind = []
        send_rank = []
        cull_inds_to_remove = []

        if n_cull > 1:  # send/recv for fixing load balance
            # CHECK FOR RANDOMNESS ISSUES AND WHICH NODES ARE USED FOR CLONES
            for r in range(size):
                # maybe remote_r should be chosen completely randomly, rather
                # than close to task of extra culled configs
                for dr in np.array(list(zip(np.array(range(1, size)), -np.array(range(1, size))))).flatten():
                    if n_cull_of_rank[r] <= max_n_cull_per_task: # not too many that need to be culled on this rank
                        break
                    # this rank has too many to cull, must receive replacement from another node
                    remote_r = (r+dr) % size
                    if n_cull_of_rank[remote_r] < max_n_cull_per_task: # send from r+dr to r
                        n_transfer = min(n_cull_of_rank[r]-max_n_cull_per_task, max_n_cull_per_task-n_cull_of_rank[remote_r])
                        recv_rank.extend([r]*n_transfer)
                        send_rank.extend([remote_r]*n_transfer)
                        local_ind = np.where(status[r, :] == 'c_t')[0][0:n_transfer]
                        recv_ind.extend(local_ind)
                        remote_ind = np.where(status[remote_r, :] == '')[0][0:n_transfer]
                        send_ind.extend(remote_ind)
                        status[r, local_ind] = 'c_s'
                        status[remote_r, remote_ind] = 'c_t_a'
                        n_cull_of_rank[r] -= n_transfer
                        n_cull_of_rank[remote_r] += n_transfer

        # save local random state, and switch to common one
        rng.switch_to_common()

        # select clones
        for r in range(size):
            list_clone_target = np.where(status[r, :] == 'c_t')[0]
            # assign clones
            n_remaining_clones = len(list_clone_target)
            while n_remaining_clones > 0:
                remote_r = rng.int_uniform(0, size)
                n_avail_remote = sum(status[remote_r, :] == '')
                if n_avail_remote > 0:  # something is available on remote_r
                    # send from random avail walker on remote_r to clone_target on r
                    n_transfer = min(n_remaining_clones, n_avail_remote)

                    # set ranks
                    send_rank.extend([remote_r]*n_transfer)
                    recv_rank.extend([r]*n_transfer)

                    # set indices
                    r_is = []
                    for ii in range(n_transfer):
                        r_i = rng.int_uniform(0, n_walkers)
                        while status[remote_r, r_i] != '':
                            r_i = rng.int_uniform(0, n_walkers)
                        # now r_i should be something with status ''
                        status[remote_r, r_i] = 'c_s'
                        r_is.append(r_i)
                    send_ind.extend(r_is)

                    status[r, list_clone_target[0:n_transfer]] = 'c_t_a'
                    recv_ind.extend(list_clone_target[0:n_transfer])

                    if n_transfer < len(list_clone_target):
                        list_clone_target = list_clone_target[n_transfer:]
                    n_remaining_clones -= n_transfer

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE POST_LOC_CLONE 15 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE POST_LOC_CLONE 16 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X POST_LOC_CLONE 17 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        # make into numpy arrays so that mathematical operations will work
        send_rank = np.array(send_rank)
        send_ind = np.array(send_ind)
        recv_rank = np.array(recv_rank)
        recv_ind = np.array(recv_ind)

        if ns_args['debug'] >= 10:
            if rank == 0:
                for i in range(len(send_rank)):
                    print(print_prefix, "send from ", send_rank[i], send_ind[i],
                          " to ", recv_rank[i], recv_ind[i])

        # save new common state, and restore to local state
        rng.switch_to_local()

        if n_cull == 1:
            if send_rank[0] == recv_rank[0] and send_rank[0] == rank:  # local copy
                walkers[recv_ind[0]].set_positions(walkers[send_ind[0]].get_positions())
                walkers[recv_ind[0]].set_cell(walkers[send_ind[0]].get_cell())
                if movement_args['do_velocities']:
                    walkers[recv_ind[0]].set_velocities(walkers[send_ind[0]].get_velocities())
                if movement_args['do_GMC']:
                    walkers[recv_ind[0]].arrays['GMC_direction'][:, :] = walkers[send_ind[0]].arrays['GMC_direction']
                if ns_args['n_extra_data'] > 0:
                    walkers[recv_ind[0]].arrays['ns_extra_data'][...] = walkers[send_ind[0]].arrays['ns_extra_data']
                if ns_args['swap_atomic_numbers']:
                    walkers[recv_ind[0]].set_atomic_numbers(walkers[send_ind[0]].get_atomic_numbers())
                    if movement_args['do_velocities']:
                        walkers[recv_ind[0]].set_masses(walkers[send_ind[0]].get_masses())
                if ns_args['track_configs']:
                    walkers[recv_ind[0]].info['config_ind'] = walkers[send_ind[0]].info['config_ind']
                    walkers[recv_ind[0]].info['from_config_ind'] = walkers[send_ind[0]].info['from_config_ind']
                    walkers[recv_ind[0]].info['config_ind_time'] = walkers[send_ind[0]].info['config_ind_time']
                walkers[recv_ind[0]].info['ns_energy'] = eval_energy(walkers[recv_ind[0]])
                if ns_args['debug'] >= 10 and size <= 1:
                    walkers[recv_ind[0]].info['n_walks'] = 0
            else:  # need send/recv
                n_send = 3*(n_atoms + 3)
                if movement_args['do_velocities']:
                    n_send += 3*n_atoms
                if movement_args['do_GMC']:
                    n_send += 3*n_atoms
                if ns_args['n_extra_data'] > 0:
                    n_send += ns_args['n_extra_data']*n_atoms
                if ns_args['swap_atomic_numbers']:
                    n_send += n_atoms  # Z
                    if movement_args['do_velocities']:
                        n_send += n_atoms  # mass
                if ns_args['track_configs']:
                    n_send += 3
                buf = np.zeros(n_send)
                if send_rank[0] == rank:  # only one config is sent/received
                    buf_o = 0
                    buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].get_positions().reshape((3*n_atoms)); buf_o += 3*n_atoms
                    buf[buf_o:buf_o+3*3] = walkers[send_ind[0]].get_cell().reshape((3*3)); buf_o += 3*3
                    if movement_args['do_velocities']:
                        buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].get_velocities().reshape((3*n_atoms)); buf_o += 3*n_atoms
                    if movement_args['do_GMC']:
                        buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].arrays['GMC_direction'].reshape((3*n_atoms)); buf_o += 3*n_atoms
                    if ns_args['n_extra_data'] > 0:
                        buf[buf_o:buf_o+ns_args['n_extra_data']*n_atoms] = walkers[send_ind[0]].arrays['ns_extra_data'].reshape((ns_args['n_extra_data']*n_atoms)); buf_o += ns_args['n_extra_data']*n_atoms
                    if ns_args['swap_atomic_numbers']:
                        buf[buf_o:buf_o+n_atoms] = walkers[send_ind[0]].get_atomic_numbers(); buf_o += n_atoms
                        if movement_args['do_velocities']:
                            buf[buf_o:buf_o+n_atoms] = walkers[send_ind[0]].get_masses(); buf_o += n_atoms
                    if ns_args['track_configs']:
                        buf[buf_o] = walkers[send_ind[0]].info['config_ind']; buf_o += 1
                        buf[buf_o] = walkers[send_ind[0]].info['from_config_ind']; buf_o += 1
                        buf[buf_o] = walkers[send_ind[0]].info['config_ind_time']; buf_o += 1
                    comm.Send([buf,  MPI.DOUBLE], dest=recv_rank[0], tag=100)
                elif recv_rank[0] == rank:
                    comm.Recv([buf, MPI.DOUBLE], source=send_rank[0], tag=100)
                    buf_o = 0
                    walkers[recv_ind[0]].set_positions(buf[buf_o:buf_o+3*n_atoms].reshape((n_atoms, 3))); buf_o += 3*n_atoms
                    walkers[recv_ind[0]].set_cell(buf[buf_o:buf_o+3*3].reshape((3, 3))); buf_o += 3*3
                    if movement_args['do_velocities']:
                        walkers[recv_ind[0]].set_velocities(buf[buf_o:buf_o+3*n_atoms].reshape((n_atoms, 3))); buf_o += 3*n_atoms
                    if movement_args['do_GMC']:
                        walkers[recv_ind[0]].arrays['GMC_direction'][:, :] = buf[buf_o:buf_o+3*n_atoms].reshape((n_atoms, 3)); buf_o += 3*n_atoms
                    if ns_args['n_extra_data'] > 0:
                        walkers[recv_ind[0]].arrays['ns_extra_data'][...] = buf[buf_o:buf_o+3*n_atoms].reshape(walkers[recv_ind[0]].arrays['ns_extra_data'].shape); buf_o += ns_args['n_extra_data']*n_atoms
                    if ns_args['swap_atomic_numbers']:
                        walkers[recv_ind[0]].set_atomic_numbers(buf[buf_o:buf_o+n_atoms].astype(int)); buf_o += n_atoms
                        if movement_args['do_velocities']:
                            walkers[recv_ind[0]].set_masses(buf[buf_o:buf_o+n_atoms]); buf_o += n_atoms
                    if ns_args['track_configs']:
                        walkers[recv_ind[0]].info['config_ind'] = int(buf[buf_o]); buf_o += 1
                        walkers[recv_ind[0]].info['from_config_ind'] = int(buf[buf_o]); buf_o += 1
                        walkers[recv_ind[0]].info['config_ind_time'] = int(buf[buf_o]); buf_o += 1
                    walkers[recv_ind[0]].info['ns_energy'] = eval_energy(walkers[recv_ind[0]])

        else:  # complicated construction of sending/receiving buffers
            # figure out how much is sent per config
            n_data_per_config = 1+3*(n_atoms + 3)
            if movement_args['do_velocities']:
                n_data_per_config += 3*n_atoms
            if movement_args['do_GMC']:
                n_data_per_config += 3*n_atoms
            if ns_args['n_extra_data'] > 0:
                n_data_per_config += ns_args['n_extra_data']*n_atoms
            if ns_args['swap_atomic_numbers']:
                n_send += n_atoms # Z
                if movement_args['do_velocities']:
                    n_send += n_atoms # mass
            if ns_args['track_configs']:
                n_data_per_config += 3

            # figure out send counts
            send_count = [0] * size
            for i in np.where(send_rank == rank)[0]:
                r_recv = recv_rank[i]
                send_count[r_recv] += n_data_per_config

            # figure out send displacements
            send_displ = [0] * size
            send_displ[0] = 0
            for i in range(1,size):
                send_displ[i] = send_displ[i-1] + send_count[i-1]

            # create empty buffer for sending
            send_count_tot = sum(send_count)
            send_data = np.zeros(send_count_tot)

            # copy data to be sent to buffer
            send_displ_t = list(send_displ)
            for i in np.where(send_rank == rank)[0]:
                r_recv = recv_rank[i]
                i_send = send_ind[i]

                data_o = send_displ_t[r_recv]
                send_data[data_o] = walkers[i_send].info['ns_energy']; data_o += 1
                send_data[data_o:data_o+3*n_atoms] = walkers[i_send].get_positions().reshape( (3*n_atoms) ); data_o += 3*n_atoms
                send_data[data_o:data_o+3*3] = walkers[i_send].get_cell().reshape( (3*3) ); data_o += 3*3
                if movement_args['do_velocities']:
                    send_data[data_o:data_o+3*n_atoms] = walkers[i_send].get_velocities().reshape( (3*n_atoms) ); data_o += 3*n_atoms
                if movement_args['do_GMC']:
                    send_data[data_o:data_o+3*n_atoms] = walkers[i_send].arrays['GMC_direction'].reshape( (3*n_atoms) ); data_o += 3*n_atoms
                if ns_args['n_extra_data'] > 0:
                    send_data[data_o:data_o+ns_args['n_extra_data']*n_atoms] = walkers[i_send].arrays['ns_extra_data'].reshape( (ns_args['n_extra_data']*n_atoms) ); data_o += ns_args['n_extra_data']*n_atoms
                if ns_args['swap_atomic_numbers']:
                    send_data[data_o:data_o+n_atoms] = walkers[i_send].get_atomic_numbers(); data_o += n_atoms
                    if movement_args['do_velocities']:
                        send_data[data_o:data_o+n_atoms] = walkers[i_send].get_masses(); data_o += n_atoms
                if ns_args['track_configs']:
                    send_data[data_o] = walkers[i_send].info['config_ind']; data_o += 1
                    send_data[data_o] = walkers[i_send].info['from_config_ind']; data_o += 1
                    send_data[data_o] = walkers[i_send].info['config_ind_time']; data_o += 1
                send_displ_t[r_recv] = data_o

            # figure out recv counts
            recv_count = [0] * size
            for i in np.where(recv_rank == rank)[0]:
                r_send = send_rank[i]
                recv_count[r_send] += n_data_per_config

            # figure out recv displacements
            recv_displ = [0] * size
            recv_displ[0] = 0
            for i in range(1,size):
                recv_displ[i] = recv_displ[i-1] + recv_count[i-1]

            # create empty buffer for receiving
            recv_count_tot = sum(recv_count)
            recv_data = np.zeros(recv_count_tot)

            # do communications
            if comm is not None:
                send_buf = [send_data, send_count, send_displ, MPI.DOUBLE]
                recv_buf = (recv_data, recv_count, recv_displ, MPI.DOUBLE)
                comm.Alltoallv(send_buf, recv_buf)
            else:
                recv_data = send_data.copy()

            # copy data from recv buffer to walkers
            recv_displ_t = list(recv_displ)
            for i in np.where(recv_rank == rank)[0]:
                r_send = send_rank[i]
                i_recv = recv_ind[i]

                data_o = recv_displ_t[r_send]
                walkers[i_recv].info['ns_energy'] = recv_data[data_o]; data_o += 1
                walkers[i_recv].set_positions( recv_data[data_o:data_o+3*n_atoms].reshape( (n_atoms, 3) )); data_o += 3*n_atoms
                walkers[i_recv].set_cell( recv_data[data_o:data_o+3*3].reshape( (3, 3) )); data_o += 3*3
                if movement_args['do_velocities']:
                    walkers[i_recv].set_velocities( recv_data[data_o:data_o+3*n_atoms].reshape( (n_atoms, 3) )); data_o += 3*n_atoms
                if movement_args['do_GMC']:
                    walkers[i_recv].arrays['GMC_direction'][:,:] = recv_data[data_o:data_o+3*n_atoms].reshape( (n_atoms, 3) ); data_o += 3*n_atoms
                if ns_args['n_extra_data'] > 0:
                    walkers[i_recv].arrays['ns_extra_data'][...] = recv_data[data_o:data_o+ns_args['n_extra_data']*n_atoms].reshape( walkers[i_recv].arrays['ns_extra_data'].shape ); data_o += ns_args['n_extra_data']*n_atoms
                if ns_args['swap_atomic_numbers']:
                    walkers[i_recv].set_atomic_numbers(recv_data[data_o:data_o+n_atoms].astype(int)); data_o += n_atoms
                    if movement_args['do_velocities']:
                        walkers[i_recv].set_masses(recv_data[data_o:data_o+n_atoms]); data_o += n_masses
                if ns_args['track_configs']:
                    walkers[i_recv].info['config_ind'] = int(recv_data[data_o]); data_o += 1
                    walkers[i_recv].info['from_config_ind'] = int(recv_data[data_o]); data_o += 1
                    walkers[i_recv].info['config_ind_time'] = int(recv_data[data_o]); data_o += 1
                recv_displ_t[r_send] = data_o

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE POST_CLONE 20 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE POST_CLONE 21 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X POST_CLONE 22 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        if ns_args['track_configs']:
            # loop over _all_ clone targets and increment cur_config_ind,
            # setting appropriate configs' new config_ind as needed
            for r in range(size):
                clone_walk_ind = np.where(status[r, :] == 'c_t_a')[0]
                for i_at in clone_walk_ind:
                    if r == rank:
                        walkers[i_at].info['from_config_ind'] = walkers[i_at].info['config_ind']
                        walkers[i_at].info['config_ind'] = cur_config_ind
                        walkers[i_at].info['config_ind_time'] = i_ns_step
                    cur_config_ind += 1
        # move cloned walkers

        if (i_ns_step == start_first_iter and movement_args['full_auto_step_sizes']):
            # set initial step sizes. Performed here since this is the first time all the arrays are in place
            conf_pre=walkers[0].copy()
            conf_pre.calc = walkers[0].calc
            move_args_pre=deepcopy(movement_args)
            walk_stats_pre=walk_single_walker(conf_pre, move_args_pre, Emax_of_step, KEmax)
            delta_step_size_setting_duration = full_auto_set_stepsizes(walkers, walk_stats_pre, movement_args, comm, Emax_of_step, KEmax, size)
            total_step_size_setting_duration += delta_step_size_setting_duration
            step_size_setting_duration += delta_step_size_setting_duration
            del(walk_stats_pre)
            del(move_args_pre)
            del(conf_pre)

        sys.stdout.flush()
        # walk clone targets
        if ns_args['debug'] >= 4:
            for i in np.where(status[rank, :] == 'c_s')[0]:
                print(print_prefix, "INFO: 30 clone source ", rank, i)
        clone_walk_ind = np.where(status[rank, :] == 'c_t_a')[0]
        for i_at in clone_walk_ind:
            if ns_args['debug'] >= 4:
                print(print_prefix, "INFO: 40 WALK clone_target ", rank, i_at)
            walk_stats = walk_single_walker(walkers[i_at], movement_args,
                                            Emax_of_step, KEmax)
            walkers[i_at].info['last_walked_iter_clone'] = i_ns_step
            # if tracking all configs, save this one that has been walked
            if track_traj_io is not None:
                walkers[i_at].info['iter'] = i_ns_step
                ase.io.write(track_traj_io, walkers[i_at], format=ns_args['config_file_format'])
            #print("WALK on rank ", rank, "at iteration ", i_ns_step, " walker ", i_at )
            if ns_args['debug'] >= 10 and size <= 1:
                walkers[i_at].info['n_walks'] += movement_args['n_model_calls']
            accumulate_stats(walk_stats_adjust, walk_stats)
            accumulate_stats(walk_stats_monitor, walk_stats)
        sys.stdout.flush()

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE POST_CLONE_WALK 25 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE POST_CLONE_WALK 26 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X POST_CLONE_WALK 27 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        # check that everything that should have been changed has, and things
        # that shouldn't have, haven't
        if ns_args['debug'] >= 10:
            final_PE_loc = [eval_energy(at, do_KE=False) for at in walkers]
            final_E_loc = [eval_energy(at) for at in walkers]
            if comm is not None:
                final_PE = np.array(comm.allgather(final_PE_loc)).flatten()
                final_E = np.array(comm.allgather(final_E_loc)).flatten()
            else:
                final_PE = final_PE_loc
                final_E = final_E_loc
            if rank == 0:
                final_status = status.flatten()
                for e in initial_unchanged:
                    if e not in final_PE:
                        print("initial_PE ", initial_PE)
                        print("final_PE ", final_PE)
                        print("initial_E ", initial_E)
                        print("final_E ", final_E)
                        print("final_status ", final_status)
                        print("WARNING: energy that should have been unchanged ", e," missing from final energies")
                for e in initial_changed:
                    if e in final_PE:
                        print("initial_PE ", initial_PE)
                        print("final_PE ", final_PE)
                        print("initial_E ", initial_E)
                        print("final_E ", final_E)
                        print("final_status ", final_status)
                        print("WARNING: energy that should have been changed ", e," still there in final energies")

        # walk extras
        if not ns_args['no_extra_walks_at_all']:
                r_i = rng.int_uniform(0, n_walkers)
                # WARNING: this may select walkers for extra walks multiple
                # times, yet never re-walk ones that were walked as clone
                # targets
                while status[rank, r_i] != '' and status[rank, r_i] != 'c_s':
                    r_i = rng.int_uniform(0, n_walkers)
                if ns_args['debug'] >= 4:
                    print(print_prefix, "INFO: 50 WALK extra ", rank, r_i)
                walk_stats = walk_single_walker(walkers[r_i], movement_args,
                                                Emax_of_step, KEmax)
                walkers[r_i].info['last_walked_iter_extra'] = i_ns_step
                # if tracking all configs, save this one that has been walked
                if track_traj_io is not None:
                    walkers[i_at].info['iter'] = i_ns_step
                    ase.io.write(track_traj_io, walkers[i_at],
                                 format=ns_args['config_file_format'])
                # print("WALK EXTRA on rank ", rank, "at iteration ", i_ns_step,
                # " walker ", r_i)
                if ns_args['debug'] >= 10 and size <= 1:
                    #walkers[r_i].info['n_walks'] += movement_args['n_steps'] # LIVIA-this gives error, does not exist
                    walkers[r_i].info['n_walks'] += movement_args['atom_traj_len']
                accumulate_stats(walk_stats_adjust, walk_stats)
                accumulate_stats(walk_stats_monitor, walk_stats)

        monitored_this_step = False
        if movement_args['monitor_step_interval'] != 0 and i_ns_step % abs(movement_args['monitor_step_interval']) == abs(movement_args['monitor_step_interval'])-1:
            adjust_step_sizes(walk_stats_monitor, movement_args, comm, monitor_only=True)
            zero_stats(walk_stats_monitor, movement_args)
            monitored_this_step = True

        if movement_args['adjust_step_interval'] != 0 and i_ns_step % abs(movement_args['adjust_step_interval']) == abs(movement_args['adjust_step_interval'])-1:

            if (not movement_args['full_auto_step_sizes']):
                adjust_step_sizes(walk_stats_adjust, movement_args, comm, do_print_rate=(not monitored_this_step))
            else:
                delta_step_size_setting_duration = full_auto_set_stepsizes(walkers, walk_stats_adjust, movement_args, comm, Emax_of_step, KEmax, size)
                total_step_size_setting_duration += delta_step_size_setting_duration
                step_size_setting_duration += delta_step_size_setting_duration
            zero_stats(walk_stats_adjust, movement_args)

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE END 30 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE END 31 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X END 32 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        if ns_args['debug'] >= 30:
            for r in range(len(status)):
                print(print_prefix, ": final status ", r, [s for s in status[r, :]])

        if (rank == 0) and ((ns_args['snapshot_interval'] > 0 and i_ns_step > 0 and i_ns_step % ns_args['snapshot_interval'] == 0) or
                            (ns_args['snapshot_seq_pairs'] and i_ns_step > 1 and i_ns_step%ns_args['snapshot_interval'] == 1) or
                            (ns_args['snapshot_time'] > 0 and time.time()-last_snapshot_time > ns_args['snapshot_time'])):
            do_snapshot = True
        else:
            do_snapshot = False
        if comm is not None:
            do_snapshot = comm.bcast(do_snapshot, root=0)
        if do_snapshot:
            save_snapshot(i_ns_step)
            last_snapshot_time = time.time()
            clean_prev_snapshot(pprev_snapshot_iter)
            pprev_snapshot_iter = prev_snapshot_iter
            prev_snapshot_iter = i_ns_step

        if ns_analyzers is not None:
            for (ns_analyzer, ns_analyzer_interval) in ns_analyzers:
                if ns_analyzer_interval > 0 and (i_ns_step+1)%ns_analyzer_interval == 0:
                    ns_analyzer.analyze(walkers, i_ns_step, "NS_loop %d" % i_ns_step)
        i_ns_step += 1
        ### END OF MAIN LOOP

    # flush remaining traj configs
    for at in traj_walker_list:
        ase.io.write(traj_io, at, parallel=False, format=ns_args['config_file_format'])
    traj_io.flush()
    traj_walker_list = []

    if ns_args['E_dump_interval'] > 0:
        if comm is not None:
            E_dump_list_all = np.array(comm.allgather(E_dump_list))
        else:
            E_dump_list_all = np.array(E_dump_list)
        if rank == 0:
            for i in range(E_dump_list_all.shape[1]):
                E_dump_io.write("step %d\n" % E_dump_list_times[i])
                if len(E_dump_list_all.shape) == 3:
                    np.savetxt(E_dump_io, E_dump_list_all[:,i,:])
                else:
                    np.savetxt(E_dump_io, E_dump_list_all[i,:])
            E_dump_io.flush()

    cur_time = time.time()
    if rank == 0:
        print( "LOOP TIME total ",cur_time-initial_time-total_step_size_setting_duration, " per iter ", (cur_time-initial_time-total_step_size_setting_duration)/(i_ns_step+1))
        print( "TIME SPENT SETTING STEP SIZES total ",total_step_size_setting_duration)

    return i_ns_step-1
```
"""

# ╔═╡ 237e5efd-32fe-4a16-b24a-34c29da83a19
md"""
## The Julian way
"""

# ╔═╡ 1cf4372f-1cac-4a09-9bf8-f76c98a3828e
md"""
```Julia
function nested_sampling_loop!(liveset::AtomWalkers, 
                                ns_params::NestedSamplingParameters, 
                                n_steps::Int64, 
                                mc_routine::MCRoutine,
                                save_strategy::DataSavingStrategy)
    df = DataFrame(iter=Int[], emax=Float64[]) # initialize empty DataFrame
    for i in 1:n_steps # main loop
        write_walker_every_n(liveset.walkers[1], i, save_strategy) # save culled walker
        iter, emax, liveset, ns_params = nested_sampling_step!(liveset, ns_params, mc_routine) # ns step
        if ns_params.fail_count >= ns_params.allowed_fail_count # reset step size
            @warn "Failed to accept MC move $(ns_params.allowed_fail_coun) times in a row. Reset step size!"
            ns_params.fail_count = 0
            ns_params.step_size = ns_params.initial_step_size
        end
        if !(iter isa typeof(missing)) # push energy into DataFrame
            push!(df, (iter, emax.val))
        end
        write_df_every_n(df, i, save_strategy) # save DataFrame
        write_ls_every_n(liveset, i, save_strategy) # save liveset
    end
    return df, liveset, ns_params
end
```
"""

# ╔═╡ 48b0f871-293d-44e2-9ffa-12481f187237
md"""## Analysis"""

# ╔═╡ 3a0c16e9-6e1e-4497-afbe-772d15105d7d
md"""
Compute the change in phase-space volume $\omega_i$:

$\omega_i = \Gamma_i - \Gamma_{i-1} = \frac{1}{N+1} \left(\frac{N}{N+1}\right)^i$
"""

# ╔═╡ 58a8ec3e-ecce-4199-912c-4418c5784cc3
omega_factors(df,n_walkers::Int64) = gamma_factors(df.iter::Vector{Int64}, n_walkers::Int64) ;

# ╔═╡ ed284d3c-9896-4271-9f50-049bae933db3
gi = omega_factors(energies,120)

# ╔═╡ ce9073cf-1c01-476b-81e6-dea76d738f9d
md"""
Shift the energies to be all positive:
"""

# ╔═╡ 1aac25a8-37cc-469d-ada4-39a027c7373c
ei = energies.emax .- minimum(energies.emax)

# ╔═╡ 312afb4b-2b97-48df-8f59-215f10678c68
md"""
Partition function:

$Z(\beta) = \sum_i \omega_i e^{-E_i \beta}$
"""

# ╔═╡ 0eed5c28-2152-4294-877a-008912486b98
partition_function(0.5, gi, ei) 

# ╔═╡ 7e267448-9bfa-4897-a117-ed75ea5e433e
kb = 8.617333262e-5 # eV/K

# ╔═╡ 601be34c-11fd-4693-9098-0eb1e4d8b0c1
ts = 10:10:10000 

# ╔═╡ 4f2a5238-b644-445d-ac33-c521d402ed41
beta = 1 ./(kb.*ts)

# ╔═╡ c58d8710-8f71-4b52-a0cf-b14760e19afb
dof = 18

# ╔═╡ 729b7408-e46e-4134-8cf6-9e7735d94840
zs = [partition_function(b, gi, ei) for b in beta]

# ╔═╡ 22169a82-0b5a-4012-8a85-64051e4ed6e9
plot(ts, zs, xlabel="Temperature (K)", ylabel="Z")

# ╔═╡ f4264836-ff91-4d89-8d59-8f2d652271f6
md"""
Comupte the internal energy

$U(\beta) = \frac{\sum_i \omega_i E_i e^{-E_i \beta}}{\sum_i \omega_i e^{-E_i \beta}}$
"""

# ╔═╡ 72329d32-e194-4090-9f41-0950adad6748
u = [internal_energy(b, gi, ei) for b in beta]

# ╔═╡ 7756264e-7e71-4d12-8a09-baeab1ee4077
plot(ts, u, xlabel="Temperature (K)", ylabel="Internal energy")

# ╔═╡ 5337631c-2150-4032-85cb-74702d5584ec
md"""
Comupte the constant-volume heat capacity:

$C_V(\beta) = \frac{dof \cdot k_B}{2} + k_B \beta^2 \left(\frac{\sum_i \omega_i E_i^2 \exp(-E_i \beta)}{Z(\beta)} - U(\beta)^2\right)$

"""

# ╔═╡ 1b7d86b4-9d84-44a6-852f-cceea593b2e9
cvs = [cv(b, gi, ei, dof) for b in beta]

# ╔═╡ 39386f04-3539-4b1f-8dec-612b081ffa46
plot(ts, cvs, xlabel="Temperature (K)",ylabel="Heat Capacity")

# ╔═╡ 6071da62-0458-4fa3-8059-f3d337aaced2
md"""
## 
![LJ6](https://github.com/yangmr04/FreeBird-demo/blob/main/LJ6.gif?raw=true)
"""

# ╔═╡ 0bfc63b6-6dc3-4ca4-85a9-7f2028dcaeb1
md"""
## How about surfaces?
"""

# ╔═╡ 6a11ff64-0ad8-4b57-a8db-06765238bfac
begin
	surf_energies = read_output("08_surf.csv")
	surf_gi = omega_factors(surf_energies, 640)
	surf_ei = surf_energies.emax .- minimum(surf_energies.emax)
	surf_ts = 1:2500
	surf_beta = 1 ./(kb.*surf_ts)
	surf_cvs = [cv(b, surf_gi, surf_ei, 3*8) for b in surf_beta]
end

# ╔═╡ 27ff4b7e-a10f-4d2c-accf-cbcd4aa2b24e
plot(surf_ts, surf_cvs, xlabel="Temperature (K)",ylabel="Heat Capacity")

# ╔═╡ 5cb4cb4e-1d02-4057-bccd-7e96e0f68d85
md"""
![surf](https://github.com/yangmr04/FreeBird-demo/blob/main/08.gif?raw=true)
"""

# ╔═╡ 7a9cb6fd-3e36-4274-8384-7abc47b22461
md"""
# Fin!

Thank you for your attention 😄
"""

# ╔═╡ Cell order:
# ╟─ad6b7485-5af8-4cd3-867a-dc3e53cf9cc0
# ╠═82211f2c-c971-48f8-a7a3-e980da7e53c0
# ╟─4e860050-5b3a-494a-8ce8-e7b5a3ebedf2
# ╟─e32b01e1-9c4d-4dd9-85e9-286fd16ed0fa
# ╟─42a7c952-e8e5-43ff-8247-520d6fad4be8
# ╟─8941699d-714b-4072-aadd-f15fb2e4adb4
# ╟─c61a733f-8a45-4481-a866-7b07d069cfd3
# ╟─3187ac88-e235-48fb-ad52-12b3fd7709ea
# ╟─ff1cf1c9-c2c9-4218-8b9d-a8c826c7d898
# ╟─e8ad89be-bbc3-4c24-800c-b3dd1400c82d
# ╟─25148d81-f849-4379-bb47-a0f26af0e257
# ╟─fe5a1609-c37e-4d75-bbe5-813454091d83
# ╠═2216fc0e-d8c3-4cfa-a003-f54701fcea66
# ╟─161a24f5-e5f1-4174-8236-ccb57d75063c
# ╟─257840d3-1004-418b-b732-200bc65fb3a0
# ╟─0d37c8c6-f858-4e58-bf41-aae0a6bee7ad
# ╟─c74d16b3-8424-4ba3-bb0c-d7434cabddcc
# ╠═9e91b9d0-6f3c-4a12-8bab-9c8b6f28857e
# ╟─de82ef79-26da-4863-bbee-82e65d01be80
# ╟─391aa43e-a4b3-4e23-a4eb-55c34af2cd51
# ╟─96978dc5-fa1e-4947-8cab-b969a67468a6
# ╠═205e7cf1-c435-4433-8f9a-9f833197b5e4
# ╟─552a77eb-f578-4497-a580-674ad9702d44
# ╟─ca3cd8bd-3057-4624-b834-35d4e76545b6
# ╠═3ff30da5-b855-4b74-9d6f-81c7f5f610c4
# ╟─9418bea4-cd5b-47eb-97fd-d3e03ced0a8c
# ╟─30095fa5-cc4b-4260-892e-f5defa326a3c
# ╠═eed6a262-f804-419d-ba96-9db596202508
# ╠═8a2920e2-9058-4e29-8560-8a249530178f
# ╠═f8d5538f-0621-4c11-875d-729af53123da
# ╠═6999a230-82d1-480d-a0fc-a7bdada6fdb6
# ╟─354c14c8-df26-460e-bb9b-c3945561d81d
# ╠═966a5497-2b72-49bf-9b76-20279a795d22
# ╠═68f4f490-f80a-4269-9000-e43ae269544a
# ╠═1beb8de4-bf20-44d3-8ff8-19a5dcfa6f94
# ╠═f2394846-280f-4c1e-9a82-5c78711452d7
# ╟─fa2a0771-d5b5-44fd-aa50-69a185377aaf
# ╟─c511a108-5b7c-4b45-8789-d6f6d55d583e
# ╟─27f32bc0-b10c-4d33-b5d6-46393e3f1a24
# ╟─237e5efd-32fe-4a16-b24a-34c29da83a19
# ╟─1cf4372f-1cac-4a09-9bf8-f76c98a3828e
# ╟─48b0f871-293d-44e2-9ffa-12481f187237
# ╟─3a0c16e9-6e1e-4497-afbe-772d15105d7d
# ╟─58a8ec3e-ecce-4199-912c-4418c5784cc3
# ╠═ed284d3c-9896-4271-9f50-049bae933db3
# ╟─ce9073cf-1c01-476b-81e6-dea76d738f9d
# ╠═1aac25a8-37cc-469d-ada4-39a027c7373c
# ╟─312afb4b-2b97-48df-8f59-215f10678c68
# ╠═0eed5c28-2152-4294-877a-008912486b98
# ╟─7e267448-9bfa-4897-a117-ed75ea5e433e
# ╟─601be34c-11fd-4693-9098-0eb1e4d8b0c1
# ╠═4f2a5238-b644-445d-ac33-c521d402ed41
# ╟─c58d8710-8f71-4b52-a0cf-b14760e19afb
# ╠═729b7408-e46e-4134-8cf6-9e7735d94840
# ╟─22169a82-0b5a-4012-8a85-64051e4ed6e9
# ╟─f4264836-ff91-4d89-8d59-8f2d652271f6
# ╠═72329d32-e194-4090-9f41-0950adad6748
# ╟─7756264e-7e71-4d12-8a09-baeab1ee4077
# ╟─5337631c-2150-4032-85cb-74702d5584ec
# ╠═1b7d86b4-9d84-44a6-852f-cceea593b2e9
# ╟─39386f04-3539-4b1f-8dec-612b081ffa46
# ╟─6071da62-0458-4fa3-8059-f3d337aaced2
# ╟─0bfc63b6-6dc3-4ca4-85a9-7f2028dcaeb1
# ╟─6a11ff64-0ad8-4b57-a8db-06765238bfac
# ╟─27ff4b7e-a10f-4d2c-accf-cbcd4aa2b24e
# ╟─5cb4cb4e-1d02-4057-bccd-7e96e0f68d85
# ╟─7a9cb6fd-3e36-4274-8384-7abc47b22461
