

function initRelatxation(system, paras)
    unitcell = [paras.L, paras.L]
    positions = copy(system.positions)

   relSys = PeriodicSystem(
        xpositions=positions,
        cutoff=system.cutoff,
        unitcell=unitcell,
        output=similar(positions),
        output_name=:forces,
    )
    return relSys
end


function update_relaxation!(x, y, i, j, d2, forces, cutoff)
    r = y - x
    dudr = 10^6 * 4 * r * (d2 - cutoff^2)

    forces[i] += dudr
    forces[j] -= dudr
    return forces
end

function update_relaxation_LJ!(x, y, i, j, d2, forces, σ,cutoff)
    r = x - y
    d = sqrt(d2)
    term2 = (σ^6 / d^6)
    term6 = (σ / d)^6
    dudr = 24  * r * (2 * term2 - term6)

    forces[i] += dudr
    forces[j] -= dudr
    return forces
end


function runRelaxation!(system, paras; nsteps=5_000)
    # initial velocities
    velocities = [ randn(eltype(system.positions)) for _ in 1:length(system.positions) ]
    dt = 1e-4
    σ = paras.sigma
    # trajectory = typeof(system.positions)[]
    for step in 1:nsteps
        # compute forces at this step
        # map_pairwise!(
        #     (x, y, i, j, d2, forces) -> update_relaxation!(x, y, i, j, d2, forces, system.cutoff),
        #     system
        # )
        map_pairwise!(
            (x, y, i, j, d2, forces) -> update_relaxation_LJ!(x, y, i, j, d2, forces, σ, system.cutoff),
            system
        )
        # Update positions and velocities
        for i in eachindex(system.positions, system.forces)
            f = system.forces[i]
            x = system.positions[i]
            v = velocities[i]
            x = x + v * dt + (f / 2) * dt^2  + sqrt(0.01* dt).*randn(2)
            v = v + f * dt
            # wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)
            # !!! IMPORTANT: Update arrays of positions and velocities
            system.positions[i] = x
            velocities[i] = v
        end
    end
end

function Relaxation!(system, paras)
    ref_pt = paras.L
    # sigma = paras.sigma
    relaxSys = initRelatxation(system, paras)
    runRelaxation!(relaxSys, paras; nsteps=5_000)
    system.positions .= relaxSys.positions
    # system.positions .= wrap_relative_to(x, SVector(ref_pt, ref_pt), system.unitcell)
    # system.positions .= copy(relaxSys.positions)
end