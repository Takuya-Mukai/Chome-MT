"""
Define interaction type

"""
mutable struct ForcesAndTorques
    force::Vector{SVector{2, Float64}}
    torque::Vector{Float64}
end

copy_output(x::ForcesAndTorques) = ForcesAndTorques(copy(x.force), copy(x.torque))

function reset_output!(output::ForcesAndTorques)
    for i in eachindex(output.force)
        output.force[i] = SVector(0.0, 0.0)
    end
    for i in eachindex(output.torque)
        output.torque[i] = 0.0
    end
    return output
end

function reducer(x::ForcesAndTorques, y::ForcesAndTorques)
    x.force .+= y.force
    x.torque .+= y.torque
    return ForcesAndTorques(x.force, x.torque)
end



"""
Function that updates forces and torque for each pair

"""
function update_interaction!(x, y, i, j, d2, orientaions, paras, output::ForcesAndTorques)
    k = paras.epsilon
    J = paras.J
    L = paras.L
    #* calculate pariwise force interaction from soft-core collision
    r = y - x
    d = sqrt(d2)
    r_hat = r ./ d
    #print("$d \n")

    # ci = (k)*x[1]/L 
    # cj = (k)*y[1]/L 
    r0 = SA[L/2, L/2]
    # rx = norm(x-r0)
    # ry = norm(y-r0)
    ci = Gaussian2D(x[1], x[2], L/2, L/2, L/4)
    cj = Gaussian2D(y[1], y[2], L / 2, L / 2, L / 4)

    c = (ci + cj) * 0.5


    # Ji = (J) * x[1] / L
    # Jj = (J) * y[1] / L

    Ji =  Gaussian2D(x[1], x[2], L / 2, L / 2, L / 4)
    Jj =  Gaussian2D(y[1], y[2], L / 2, L / 2, L / 4)


    Jij = (Ji+Jj) * 0.5
    # c = 1.0
    # f_ij = 2c * r ./ d
    σ = paras.sigma
    # F_LJ = SA[0., 0.]
    # if d < σ * 0.05 
        # F_LJ =  -(σ - d)^2 * r ./d
        # F_LJ = SA[0., 0.]
    # else
        # F_LJ = c/d2 * r 
    # end

    
    # σ = paras.sigma
    # term2 = (σ^6 / d^6)
    # term6 = (σ / d)^6
    # F_LJ = 24 *c * r* (2*term2 - term6)

    if d < 2σ
        F_ve = -1*(2σ - d) * r_hat
        F_at = SA[0.0, 0.0]
    else
        F_ve = SA[0.0, 0.0]
        # F_at = c * r_hat  / d2
        term6 = (2σ / d)^6
        F_at = 24 * c * r * (- term6)
    end
    
    
    # if d == 0.0
    #     @show x, y, F_LJ
    # end
    # end
    # @show x,y
    # if d < 1e-3
    #     @show x, y 
    #     F_LJ = SA[0.,0.]

    output.force[i] += (F_at .+ F_ve)
    output.force[j] -= (F_at .+ F_ve)
    
    # output.force[i] += (F_LJ) 
    # output.force[j] -= (F_LJ)

    #* calculate torque for each pari n_i x n_j
    # if d < cutoff
    n_cross =   -Jij*sin(orientaions[i] - orientaions[j]) / d
    # else
    #     n_cross = SVector(0.0, 0.0, 0.0)
    #     print("this should not be printed if the package cuttoff the interaction")
    # end

    output.torque[i] += n_cross
    output.torque[j] -= n_cross

    return output
end
