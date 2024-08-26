"""
Function that updates pariwise forces and torque for each pair in spatio dependance Viscek type interaction

"""
function update_interaction!(x, y, i, j, d2, orientaions, paras, spatio_func, output::ForcesAndTorques)
    # println("here")
    # ϵ = paras.epsilon
    J = paras.J
    # L = paras.L
    σ = paras.sigma
    
    Jij = J*spatio_func(x[2],y[2], paras.L)
    #* calculate pariwise force interaction from hardcore repulsion
    # r = y - x
    d = sqrt(d2)
    # r_hat = r ./ d
    Fij = SA[0.0, 0.0]
    # if (d < 2σ ) < 0.0
    #     Fij = -ϵ*(2σ-d)*r_hat
    # end

    output.force[i] += (Fij)
    output.force[j] -= (Fij)

    #* calculate pariwise alignment
    Tij = 0.0
    if d < 1.0
        Tij = - Jij * sin(orientaions[i] - orientaions[j])
    end
    output.torque[i] += Tij
    output.torque[j] -= Tij

    return output
end
