


function initPosition(paras, N)
    positions = [SA[0., 0.] for _ in 1:N]
    nL = sqrt(N)
    dx = paras.L / sqrt(N)
    for i in 1:N
        x = dx*(i % nL) +randn()*paras.sigma*0.3
        y = floor(i / nL) * dx + randn()*paras.sigma*0.1
        positions[i] = SA[x,y]
    end
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


function simulate(system, paras; nsteps::Int=1000, dt::Float64=2e-4)

    # Relaxation!(system, paras)

    N = length(system.positions)
    L = paras.L
    Dr = paras.Dr
    ref_pos = L / 2
    velocities = [ randn(eltype(system.positions)) for _ in 1:length(system.positions) ]
    orientaions = rand(N) * 2π
    vel = paras.vel

    all_pos = [similar(system.positions) for _ in 1:nsteps]
    all_ϕ = [similar(orientaions) for _ in 1:nsteps]
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

        all_pos[step] = copy(system.positions)
        all_ϕ[step] = copy(orientaions)
        # end
    end
    return all_pos, all_ϕ
end


function normDistri(x; sigma=1, mu = 0)
    return  1/(sigma*sqrt(2π)) * exp(-0.5*((x-mu)/sigma)^2)
end



function Gaussian2D(x, y, x0, y0, σ)
    exp(-((x-x0)^2 / (2σ^2) + (y-y0)^2 / (2σ^2)))
end