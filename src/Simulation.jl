using Base.Threads


function initPosition(paras, N)
    positions = [SA[0., 0.] for _ in 1:N]
    nL = sqrt(N)
    dx = paras.L / sqrt(N)
    # for i in 1:N
    #     x = dx*(i % nL) +randn()*paras.sigma*0.3
    #     y = floor(i / nL) * dx + randn()*paras.sigma*0.1
    #     positions[i] = SA[x,y]
    # end

    positions = [SA[rand(), rand()] .- 0.5 for _ in 1:N] .* paras.L .+ (SA[paras.L/2, paras.L/2], )
    # positions = [SA[paras.L/30, paras.L/30] .+ randn(2)]
    return positions
end

function initSystem(paras::Paras; N::Int=10)
    # if isnothing(positions)
    #     positions = randUnitVec(3, N) .* paras.R
    # end
    unitcell = [paras.L, paras.L]
    # positions = [SA[rand()*unitcell[1], rand()*unitcell[2]] for _ in 1:N]
    positions = initPosition(paras, N)
    system = PeriodicSystem(
        xpositions=positions,
        cutoff=paras.cutoff,
        unitcell=unitcell,
        output=ForcesAndTorques(similar(positions), @SVector zeros(N)),
        output_name=:ForcesAndTorques
    )

    return system
end


function simulate(system, paras; nsteps::Int=1000, dt::Float64=2e-4, isave = 500)

    # Relaxation!(system, paras)

    N = length(system.positions)
    L = paras.L
    Dr = paras.Dr
    ref_pos = L / 2
    velocities = [ randn(eltype(system.positions)) for _ in 1:length(system.positions) ]
    orientaions = rand(N) * 2π
    vel = paras.vel

    # all_pos = [similar(system.positions) for _ in 1:nsteps]
    # all_ϕ = [similar(orientaions) for _ in 1:nsteps]

    all_pos = typeof(system.positions)[]
    all_ϕ = typeof(orientaions)[]

    # force_cache = [similar(system.F)]
    for step in 1:nsteps
        # if findall(x->x==SA[0.,0.], positions) == 0
        #     prinln("help!")
        # end
        # compute forces at this step
        map_pairwise!(
            (x, y, i, j, d2, output) -> update_interaction!(x, y, i, j, d2, orientaions, paras, output),
            system
        )


        #* Update positions and velocities
        @inbounds  for i in eachindex(system.positions)
            f = system.ForcesAndTorques.force[i]
            x = system.xpositions[i]
            ϕ = orientaions[i]
            # v = velocities[i]
            # vel = 0.8*x[1] / L + 0.2
            x = x + SA[vel*cos(ϕ), vel*sin(ϕ)] * dt + f * dt

            # x = x + v * dt + (f / 2) * dt^2 + vel * SA[cos(ϕ), sin(ϕ)] * dt
            # v = v + f * dt

            ϕ = ϕ + system.ForcesAndTorques.torque[i] * dt + sqrt(2 * Dr * dt) * randn()

            #* wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(ref_pos, ref_pos), system.unitcell)
            # x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)

            # !!! IMPORTANT: Update arrays of positions and velocities
            system.positions[i] = x
            orientaions[i] = ϕ
            # velocities[i] = v
        end

        # all_pos[step] = copy(system.positions)
        # all_ϕ[step] = copy(orientaions)
        if step % isave == 0
            push!(all_pos, copy(system.positions))
            push!(all_ϕ, copy(orientaions))
        end
        # end
    end
    return all_pos, all_ϕ
end

function simulate_midStep(system, paras; nsteps::Int=1000, dt::Float64=2e-4, isave=500)

    # Relaxation!(system, paras)

    N = length(system.positions)
    L = paras.L
    Dr = paras.Dr
    ref_pos = L / 2
    velocities = [randn(eltype(system.positions)) for _ in 1:length(system.positions)]
    orientaions = rand(N) * 2π
    vel = paras.vel

    # all_pos = [similar(system.positions) for _ in 1:nsteps]
    # all_ϕ = [similar(orientaions) for _ in 1:nsteps]

    all_pos = typeof(system.positions)[]
    all_ϕ = typeof(orientaions)[]

    # force_cache = [similar(system.F)]

    for step in 1:nsteps
        posisitons_tmp = copy(system.positions)
        orienations_tmp = copy(orientaions)
        # if findall(x->x==SA[0.,0.], positions) == 0
        #     prinln("help!")
        # end
        # compute forces at this step

        # !!! mid step
        map_pairwise!(
            (x, y, i, j, d2, output) -> update_interaction_FEND!(x, y, i, j, d2, orientaions, paras, output),
            system
        )

        #* Update positions and velocities in mid step
        # eta = @SVector randn(N)
        @inbounds for i in eachindex(system.positions)
            f_mid = system.ForcesAndTorques.force[i]
            x_mid = system.xpositions[i]
            ϕ_mid = orientaions[i]
            # v = velocities[i]
            # vel = 0.8*x[1] / L + 0.2
            x_mid = x_mid + SA[vel*cos(ϕ_mid), vel*sin(ϕ_mid)] * 0.5*dt + 0.5*f_mid * dt

            # x = x + v * dt + (f / 2) * dt^2 + vel * SA[cos(ϕ), sin(ϕ)] * dt
            # v = v + f * dt

            ϕ_mid = ϕ_mid + system.ForcesAndTorques.torque[i] * 0.5*dt

            #* wrapping to origin for obtaining a pretty animation
            x_mid = wrap_relative_to(x_mid, SVector(ref_pos, ref_pos), system.unitcell)
            # x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)

            # !!! IMPORTANT: Update arrays of positions and velocities
            system.positions[i] = x_mid
            orientaions[i] = ϕ_mid
            # velocities[i] = v
        end

        # !!! final step
        map_pairwise!(
            (x, y, i, j, d2, output) -> update_interaction_FEND!(x, y, i, j, d2, orientaions, paras, output),
            system
        )
        
        @inbounds  for i in eachindex(system.positions)
            f = system.ForcesAndTorques.force[i]
            # x = system.xpositions[i]
            x = posisitons_tmp[i]
            ϕ = orienations_tmp[i]
            # v = velocities[i]
            # vel = 0.8*x[1] / L + 0.2
            x = x + SA[vel*cos(ϕ), vel*sin(ϕ)] * dt + f * dt

            # x = x + v * dt + (f / 2) * dt^2 + vel * SA[cos(ϕ), sin(ϕ)] * dt
            # v = v + f * dt

            ϕ = ϕ + system.ForcesAndTorques.torque[i] * dt + sqrt(2 * Dr * dt) * randn()

            #* wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(ref_pos, ref_pos), system.unitcell)
            # x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)

            # !!! IMPORTANT: Update arrays of positions and velocities
            # velocities[i] = v
            system.positions[i] = x
            orientaions[i] = ϕ
        end

        # all_pos[step] = copy(system.positions)
        # all_ϕ[step] = copy(orientaions)
        if step % isave == 0
            push!(all_pos, copy(system.positions))
            push!(all_ϕ, copy(orientaions))
        end
        # end
    end
    return all_pos, all_ϕ
end



function Gaussian1D(x; sigma=1, mu=0)
    return  1/(sigma*sqrt(2π)) * exp(-0.5*((x-mu)/sigma)^2)
end



function Gaussian2D(x, y, x0, y0, σ)
    exp(-((x-x0)^2 / (2σ^2) + (y-y0)^2 / (2σ^2)))
end



function simulation_chemotaxis(system, paras; nsteps::Int=1000, dt::Float64=2e-4, isave = 100)
    N = length(system.positions)
    L = paras.L 
    Dr = paras.Dr
    ref_pos = L/2
    dx = paras.dx
    ω0 = 0.0
    ω = 3.
    α = 7.0

    velocities = [randn(eltype(system.positions)) for _ in 1:length(system.positions)]
    orientaions = rand(N) * 2π
    vel = paras.vel

    all_pos = typeof(system.positions)[]
    all_ϕ = typeof(orientaions)[]

    # all_field = typeof(field)[]

    # chemo_force = similar(system.ForcesAndTorques.force)
    # chemo_torque = similar(system.ForcesAndTorques.torque)


    nx = round(Int, L / dx) 
    @show nx

    field = zeros(nx, nx)
    # field = [(x * L)/nx for x in 1:nx, y in 1:nx]
    field0 = copy(field)

    all_field = typeof(field)[]

    colloild_bound = BoundVec(paras.sigma, dx)


    for step in 1:nsteps
        # calculate pariwise force
        map_pairwise!(
            (x, y, i, j, d2, output) -> update_interaction!(x, y, i, j, d2, orientaions, paras, output),
            system
        )

        """updata position and orientaions"""
        @inbounds for i in eachindex(system.xpositions)
            x = system.xpositions[i]
            ϕ = orientaions[i]

            F_chem = chemotaxis_periodic(x, field0, paras, colloild_bound)

            f = system.ForcesAndTorques.force[i] - α * F_chem

            x = x + SA[vel*cos(ϕ), vel*sin(ϕ)] * dt + f * dt + sqrt(2 * Dr * dt) .* randn(2)

            T = ω0 +  system.ForcesAndTorques.torque[i] + ω * cross(SA[cos(ϕ), sin(ϕ)], F_chem) 
            # T = ω0 + system.ForcesAndTorques.torque[i] + ω * sin(atan(F_chem[2], F_chem[1]) - ϕ) + sqrt(2 * Dr * dt) * randn()
            ϕ = ϕ + T * dt + sqrt(2 * Dr * dt) * randn()

            # !!! IMPORTANT: Update arrays of positions and velocities
            x = wrap_relative_to(x, SVector(ref_pos, ref_pos), system.unitcell)
            system.positions[i] = x
            orientaions[i] = ϕ
        end

        """
        updata field
        """
        field, field0 = UpdateField(system.xpositions, field, field0, paras, dt)
        field, field0 = field0, field

        #* save data
        if step % isave == 0
            push!(all_pos, copy(system.positions))
            push!(all_ϕ, copy(orientaions))
            push!(all_field, copy(field0))
        end
    end

    return all_pos, all_ϕ, all_field
end



function simulation_visceks(system, paras; nsteps::Int=1000, dt::Float64=2e-4, isave=100)
    N = length(system.positions)
    L = paras.L
    Dr = paras.Dr
    ref_pos = L / 2
    dx = paras.dx
    ω0 = 0.0
    # ω = 3.0
    # α = 7.0

    velocities = [randn(eltype(system.positions)) for _ in 1:length(system.positions)]
    orientaions = rand(N) * 2π
    vel = paras.vel

    all_pos = typeof(system.positions)[]
    all_ϕ = typeof(orientaions)[]

    # all_field = typeof(field)[]

    # chemo_force = similar(system.ForcesAndTorques.force)
    # chemo_torque = similar(system.ForcesAndTorques.torque)


    # nx = round(Int, L / dx)
    # @show nx

    # field = zeros(nx, nx)
    # field = [(x * L)/nx for x in 1:nx, y in 1:nx]
    # field0 = copy(field)

    # all_field = typeof(field)[]

    # colloild_bound = BoundVec(paras.sigma, dx)


    for step in 1:nsteps
        # calculate pariwise force
        map_pairwise!(
            (x, y, i, j, d2, output) -> update_interaction!(x, y, i, j, d2, orientaions, paras, output),
            system
        )

        """updata position and orientaions"""
        @inbounds for i in eachindex(system.xpositions)
            x = system.xpositions[i]
            ϕ = orientaions[i]

            # F_chem = chemotaxis_periodic(x, field0, paras, colloild_bound)

            f = system.ForcesAndTorques.force[i]

            x = x + SA[vel*cos(ϕ), vel*sin(ϕ)] * dt + f * dt

            T = ω0 + system.ForcesAndTorques.torque[i] 
            # T = ω0 + system.ForcesAndTorques.torque[i] + ω * sin(atan(F_chem[2], F_chem[1]) - ϕ) + sqrt(2 * Dr * dt) * randn()
            ϕ = ϕ + T * dt + sqrt(2 * Dr * dt) * randn()

            # !!! IMPORTANT: Update arrays of positions and velocities
            x = wrap_relative_to(x, SVector(ref_pos, ref_pos), system.unitcell)
            system.positions[i] = x
            orientaions[i] = ϕ
        end

        # """
        # updata field
        # """
        # field, field0 = UpdateField(system.xpositions, field, field0, paras, dt)
        # field, field0 = field0, field

        #* save data
        if step % isave == 0
            push!(all_pos, copy(system.positions))
            push!(all_ϕ, copy(orientaions))
        end
    end

    return all_pos, all_ϕ
end


function simulation_visceks_spatio(system, paras, spatio_func; nsteps::Int=1000, dt::Float64=2e-4, isave=100)
    N = length(system.positions)
    L = paras.L
    Dr = paras.Dr
    ref_pos = L / 2
    dx = paras.dx
    ω0 = 0.0
    # ω = 3.0
    # α = 7.0

    velocities = [randn(eltype(system.positions)) for _ in 1:length(system.positions)]
    orientaions = rand(N) * 2π
    vel = paras.vel

    all_pos = typeof(system.positions)[]
    all_ϕ = typeof(orientaions)[]

    # all_field = typeof(field)[]

    # chemo_force = similar(system.ForcesAndTorques.force)
    # chemo_torque = similar(system.ForcesAndTorques.torque)


    # nx = round(Int, L / dx)
    # @show nx

    # field = zeros(nx, nx)
    # field = [(x * L)/nx for x in 1:nx, y in 1:nx]
    # field0 = copy(field)

    # all_field = typeof(field)[]

    # colloild_bound = BoundVec(paras.sigma, dx)


    for step in 1:nsteps
        # calculate pariwise force
        map_pairwise!(
            (x, y, i, j, d2, output) -> update_interaction!(x, y, i, j, d2, orientaions, paras, spatio_func, output),
            system
        )

        """updata position and orientaions"""
        for i in eachindex(system.xpositions)
            x = system.xpositions[i]
            ϕ = orientaions[i]

            # F_chem = chemotaxis_periodic(x, field0, paras, colloild_bound)

            f = system.ForcesAndTorques.force[i]

            x = x + SA[vel*cos(ϕ), vel*sin(ϕ)] * dt + f * dt

            T = ω0 + system.ForcesAndTorques.torque[i] 
            # T = ω0 + system.ForcesAndTorques.torque[i] + ω * sin(atan(F_chem[2], F_chem[1]) - ϕ) + sqrt(2 * Dr * dt) * randn()
            ϕ = ϕ + T * dt + sqrt(2 * Dr * dt) * randn()

            # !!! IMPORTANT: Update arrays of positions and velocities
            x = wrap_relative_to(x, SVector(ref_pos, ref_pos), system.unitcell)
            system.positions[i] = x
            orientaions[i] = ϕ
        end

        # """
        # updata field
        # """
        # field, field0 = UpdateField(system.xpositions, field, field0, paras, dt)
        # field, field0 = field0, field

        #* save data
        if step % isave == 0
            push!(all_pos, copy(system.positions))
            push!(all_ϕ, copy(orientaions))
        end
    end

    return all_pos, all_ϕ
end
