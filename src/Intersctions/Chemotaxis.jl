"""
Function that updates pariwise forces and torque for each pair

"""
function update_interaction!(x, y, i, j, d2, orientaions, paras, output::ForcesAndTorques)
    ϵ = paras.epsilon
    J = paras.J
    L = paras.L
    σ = paras.sigma
    
    #* calculate pariwise force interaction from hardcore repulsion
    #* overlapping particles are shifted in the direction of their distance vector by equal amounts until a distance of 2σ
    r = y - x
    d = sqrt(d2)
    r_hat = r ./ d
    Fij = SA[0.0, 0.0]
    # if (d < 2σ ) < 0.0
        Fij = -ϵ*(2σ-d)*r_hat
    # end

    output.force[i] += (Fij)
    output.force[j] -= (Fij)

    #* calculate pariwise alignment
    Tij = 0.0
    if d < 1.5
        Tij = - J * sin(orientaions[i] - orientaions[j]) / d
    end
    output.torque[i] += Tij
    output.torque[j] -= Tij

    return output
end



"""
Chemotaxis interaction
"""

function chemotaxis_periodic(pos, field, paras, colloild_bound)
    σ = paras.sigma
    dx = paras.dx
    r0 = pos
    nx, ny = size(field)
    unit_vec, ϕ, dϕ = colloild_bound

    force = SA[0.0, 0.0]
    # torque = 0.0

    lo = r0 .- 1.5σ
    up = r0 .+ 1.5σ
    box_range = round(Int, 3σ/dx)
    refpoint = Array{Float64}(undef, box_range, box_range)


    # #* mapping physical position to grid index
    lo_idx, up_idx = pos2idx(lo, dx), pos2idx(up, dx)
    xlo, ylo = lo_idx
    xup, yup = up_idx
    
    xrange = lo_idx[1] : lo_idx[1] + box_range - 1
    yrange = lo_idx[2] : lo_idx[2] + box_range - 1

    xrange = wrap_index.(xrange, nx)
    yrange = wrap_index.(yrange, ny)
    for (j, yj) in enumerate(yrange)
        for (i, xi) in enumerate(xrange)
            refpoint[i,j] = field[xi, yj]
        end
    end

    xs = 0 : dx : box_range*dx-dx
    ys = 0 : dx : box_range*dx-dx
    # refpoint = [field[i,j] for i in lo_idx[1]:up_idx[1], j in lo_idx[2]:up_idx[2]]
    # xlist = round((xlo - 1) * dx, digits=3):dx:round((xup - 1) * dx, digits=3)
    # ylist = round((ylo - 1) * dx, digits=3):dx:round((yup - 1) * dx, digits=3)


    # #* interpolation
    # sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xs, ys)
    x, y = r0
    #* find the gradient across a finte size particle
    for j in 1:length(ϕ)
        # unit_vec[i, j] = SA[sin(θ[i])cos(ϕ[j]), sin(θ[i])sin(ϕ[j]), cos(θ[i])]
        n = unit_vec[j]
        pt1 = n .* σ
        # pt1 = wrap_position.(r0 + pt1, L)
        pt_in_ref = r0 + pt1 - idx2pos.(lo_idx, dx)
        # ∇C = Grad(pt, dx, field)

        ∇C = Interpolations.gradient(sitp, pt_in_ref[1], pt_in_ref[2])
        ∇C = ∇C .- dot(∇C, n) .* n

        force += ∇C * σ * dϕ
        # torque += cross(∇C, n) * σ * dϕ
    end
    A = 2π*σ 
    return -force / A
end


function chemotaxis(pos, field, paras, colloild_bound)
    σ = paras.sigma
    dx = paras.dx
    r0 = pos
    L = paras.L

    unit_vec, ϕ, dϕ = colloild_bound

    force = SA[0.0, 0.0]
    torque = 0.0

    lo = r0 .- 1.5σ
    up = r0 .+ 1.5σ

    #* mapping physical position to grid index
    lo_idx, up_idx = pos2idx(lo, dx), pos2idx(up, dx)
    xlo, ylo = lo_idx
    xup, yup = up_idx

    refpoint = [field[i, j] for i in lo_idx[1]:up_idx[1], j in lo_idx[2]:up_idx[2]]
    

    xlist = round((xlo - 1) * dx, digits=3):dx:round((xup - 1) * dx, digits=3)
    ylist = round((ylo - 1) * dx, digits=3):dx:round((yup - 1) * dx, digits=3)

    # xs = 0:dx:L-dx

    #* interpolation
    sitp = scale(interpolate(refpoint, BSpline(Quadratic(InPlace(OnCell())))), xlist, ylist)
    # sitp = scale(interpolate(field, BSpline(Quadratic(InPlace(OnCell())))), xs, xs)

    x, y = r0
    #* find the gradient across a finte size particle
    for j in 1:length(ϕ)
        # unit_vec[i, j] = SA[sin(θ[i])cos(ϕ[j]), sin(θ[i])sin(ϕ[j]), cos(θ[i])]
        n = unit_vec[j]
        pt1 = n .* σ
        ∇C = Interpolations.gradient(sitp, x + pt1[1], y + pt1[2])
        # ∇C = Interpolations.gradient(sitp, wrap_position(x + pt1[1], L), wrap_position(y + pt1[2], L))
        ∇C = ∇C .- dot(∇C, n) .* n

        force += ∇C * σ * dϕ
        torque += cross(∇C, n)
    end
    A = 2π * σ
    return -force / A
end


function BoundVec(R, dx)
    n = ceil(Int, 2π*R/dx)
    ϕ = LinRange(0, 2π - (2π / n), n)
    dϕ = ϕ[2] - ϕ[1]
    unit_vec = Array{SVector{2,Float64},1}(undef, length(ϕ))

    for j in eachindex(ϕ)
        unit_vec[j] = SA[cos(ϕ[j]), sin(ϕ[j])]
    end


    return unit_vec, ϕ, dϕ
end

# function Disk(; n=90, R=15, pos=SA[0.0, 0.0])
#     edges = [pos .+ R * SA[cos(θ), sin(θ)] for θ in LinRange(0, 2π - (2π / n), n)]
#     slip = [@SVector zeros(n) for _ in 1:2]
#     dl = R * 2π / n

#     return Disk(n, R, pos, edges, slip, dl)
# end


function Grad(pt, dx, field)
    lo = floor.(Int, pt / dx) .+ 1
    up = ceil.(Int, pt / dx) .+ 1
    nx, ny = size(field)
    lo, up = wrap_index.(lo, nx), wrap_index.(up, ny)
    x1, y1 = lo
    x2, y2 = up

    ∇x = (field[x2] - field[x1]) / dx 
    ∇y = (field[y2] - field[y1]) / dx
    
    return SA[∇x, ∇y]
end

