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



# """
# Function that updates forces and torque for each pair

# """
# function update_interaction!(x, y, i, j, d2, orientaions, paras, output::ForcesAndTorques)
#     ϵ = paras.epsilon
#     J = paras.J
#     L = paras.L
#     #* calculate pariwise force interaction from soft-core collision
#     r = y - x
#     d = sqrt(d2)
#     r_hat = r ./ d
#     #print("$d \n")

#     # ci = (k)*x[1]/L 
#     # cj = (k)*y[1]/L 
#     # r0 = SA[L/2, L/2]
#     # rx = norm(x-r0)
#     # # ry = norm(y-r0)
#     # ci = Gaussian2D(x[1], x[2], L/2, L/2, L/4)
#     # cj = Gaussian2D(y[1], y[2], L / 2, L / 2, L / 4)

#     ci = Gaussian1D(x[2], sigma=L/4, mu=L/2)
#     cj = Gaussian1D(y[2], sigma=L / 4, mu=L / 2)

#     # ci = 1.6 - Gaussian1D(x[2], sigma=L/4, mu=L/2) * 100
#     # cj = 1.6 - Gaussian1D(y[2], sigma=L / 4, mu=L / 2) * 100


#     # ci = Gaussian2D(x[1], x[2], L / 4, L / 4, L / 8) 
#     #     + Gaussian2D(x[1], x[2], L / 4, 3*L / 4, L / 8)
#     #     + Gaussian2D(x[1], x[2], 3*L / 4, L / 4, L / 8)
#     #     + Gaussian2D(x[1], x[2], 3*L / 4, 3*L / 4, L / 8)

#     # cj = Gaussian2D(x[1], x[2], L / 4, L / 4, L / 4)
#     #     +Gaussian2D(x[1], x[2], L / 4, 3 * L / 4, L / 4)
#     #     +Gaussian2D(x[1], x[2], 3 * L / 4, L / 4, L / 4)
#     #     +Gaussian2D(x[1], x[2], 3 * L / 4, 3 * L / 4, L / 4)

#     c = (ci + cj) * 0.5 * ϵ
#     # c = paras.epsilon   


#     # Ji = (J) * x[1] / L
#     # Jj = (J) * y[1] / L

#     # Ji =  Gaussian2D(x[1], x[2], L / 2, L / 2, L / 4)
#     # Jj =  Gaussian2D(y[1], y[2], L / 2, L / 2, L / 4)


#     # Jij = (Ji+Jj) * 0.5
#     Jij = J
#     # c = 1.0
#     # f_ij = 2c * r ./ d
#     σ = paras.sigma
#     # F_LJ = SA[0., 0.]
#     # if d < σ * 0.05 
#         # F_LJ =  -(σ - d)^2 * r ./d
#         # F_LJ = SA[0., 0.]
#     # else
#         # F_LJ = c/d2 * r 
#     # end

    
#     # σ = paras.sigma
#     # term2 = (σ^6 / d^6)
#     # term6 = (σ / d)^6
#     # F_LJ = 24 *c * r* (2*term2 - term6)

#     # if d < 2σ
#     #     F_ve = -1*(2σ - d) * r_hat
#     #     F_at = SA[0.0, 0.0]
#     # else
#         # F_ve = SA[0.0, 0.0]
#     #     # F_at = c * r_hat  / d2
#     #     term6 = (2σ / d)^6
#     #     F_at = 24 * c * r * (- term6)
#     # end
    
    
#     # if d == 0.0
#     #     @show x, y, F_LJ
#     # end
#     # end
#     # @show x,y
#     # if d < 1e-3
#     #     @show x, y 
#     #     F_LJ = SA[0.,0.]


#     #* log-bond 
#     fij = 0.0
#     # k = 1.0
#     if d < (2σ)
#     # #     fwca = 24 * ϵ * d2inv * (2 * (σ6 * d6inv)^2 - σ6 * d6inv)
#         fij = -1.0*(2σ-d)
#     else
#     # end
    
#         fij = - c / (1 - d2 / (2σ)^2)
#     end
#     #* string force 
#     # fij = 0.5*c*(d-2d)^2
#     # fij = -c / (1 - d2 / (2σ)^2)
#     f_fend = fij * r_hat
#     # F_at = c/d * r
#     # output.force[i] += (F_at .+ F_ve)
#     # output.force[j] -= (F_at .+ F_ve)
    
#     output.force[i] += (f_fend)
#     output.force[j] -= (f_fend)

#     #* calculate torque for each pari n_i x n_j
#     if d < (2σ*1.2)
#         n_cross =   -Jij*sin(orientaions[i] - orientaions[j])
#     else
#         n_cross = 0.0
#     end
#     #     n_cross = SVector(0.0, 0.0, 0.0)
#     #     print("this should not be printed if the package cuttoff the interaction")
#     # end

#     output.torque[i] += n_cross
#     output.torque[j] -= n_cross

#     return output
# end

# function update_interaction_FEND!(x, y, i, j, d2, orientaions, paras, output::ForcesAndTorques)
#     k = paras.epsilon
#     σ = paras.sigma
#     # J = paras.J
#     L = paras.L 
#     #* calculate pariwise force interaction from soft-core collision
#     r = y - x
#     d = sqrt(d2)
#     # r_hat = r ./ d

#     # ab = vector(coord_i, coord_j, boundary)
#     # ab = y-x
#     # r = norm(ab)
#     # d2 = d^2

#     ci = Gaussian1D(x[2], sigma=L / 4, mu=L / 2)
#     cj = Gaussian1D(y[2], sigma=L / 4, mu=L / 2)

#     c = (ci + cj) * 0.5 * k * 100
#     # c = 0.01

#     d2inv = inv(d2)
#     d6inv = d2inv^3
#     σ6 = (2σ)^6
#     fwca_divr = 0.0
#     fmag_divr = 0.0

#     if d < ((2σ) * 2^(1 / 6))
#         fwca_divr = 24 * c * d2inv * (2 * (σ6 * d6inv)^2 - σ6 * d6inv)
#         # fwca_divr = 0.0
#     end
#     fmag_divr = fwca_divr - c / (1 - d2 / (2σ*1.05)^2)
#     # fmag_divr = fwca_divr - 0.0

#     f = fmag_divr * r
#     # fend = 0.0 * r
#     # @show f
#     output.force[i] += (f)
#     output.force[j] -= (f)
#     # @show i, j output.force[i], output.force[j]
#     n_cross = 0.0
#     if d < (2σ * 1.05)
#         n_cross = -c * sin(orientaions[i] - orientaions[j])
#     else
#         n_cross = 0.0
#     end

#     output.torque[i] += n_cross
#     output.torque[j] -= n_cross

#     return output
# end

# function update_interaction_FEND!(x, y, i, j, d2, orientaions, paras, output::ForcesAndTorques)
#     k = paras.epsilon
#     σ = paras.sigma
#     J = paras.J
#     L = paras.L
#     #* calculate pariwise force interaction from soft-core collision
#     r = y - x
#     d = sqrt(d2)
#     r_hat = r ./ d

#     # ab = vector(coord_i, coord_j, boundary)
#     # ab = y-x
#     # r = norm(ab)
#     # d2 = d^2

#     # ci = Gaussian1D(x[2], sigma=L / 4, mu=L / 2)
#     # cj = Gaussian1D(y[2], sigma=L / 4, mu=L / 2)

#     # ci = Gaussian1D(x[2], sigma=L/4, mu=L/2) *100 + .05
#     # cj = Gaussian1D(y[2], sigma=L / 4, mu=L / 2)*100 + .05

#     # ci = .517 - Gaussian1D(x[2], sigma=L / 4, mu=L / 2) * 100
#     # cj = .517 - Gaussian1D(y[2], sigma=L / 4, mu=L / 2) * 100


#     ci = 1.05 - (Gaussian2D(x[1], x[2], L / 4, L / 4, L / 6) 
#         + Gaussian2D(x[1], x[2], L / 4, 3*L / 4, L / 6)
#         + Gaussian2D(x[1], x[2], 3*L / 4, L / 4, L / 6)
#         + Gaussian2D(x[1], x[2], 3*L / 4, 3*L / 4, L / 6))


#     cj = 1.05 - (Gaussian2D(x[1], x[2], L / 4, L / 4, L / 6)
#                  + Gaussian2D(x[1], x[2], L / 4, 3 * L / 4, L / 6)
#                  + Gaussian2D(x[1], x[2], 3 * L / 4, L / 4, L / 6)
#                  + Gaussian2D(x[1], x[2], 3 * L / 4, 3 * L / 4, L / 6))

#     # ci = 1.1 - Gaussian2D(x[1], x[2], L/2, L/2, L/4)
#     # cj = 1.1 - Gaussian2D(y[1], y[2], L / 2, L / 2, L / 4)

#     c = (ci + cj) * 0.5 * k
#     Jij = J * c / k

#     # ab = vector(coord_i, coord_j, boundary)
#     f_bond = 0.0
#     if d < (2σ * 1.05)
#         f_bond = c * (d - 2σ)
#     end
#     f = f_bond * r_hat

#     # f = fmag_divr * r
#     # fend = 0.0 * r
#     # @show f
#     output.force[i] += (f)
#     output.force[j] -= (f)
#     # @show i, j output.force[i], output.force[j]
#     n_cross = 0.0
#     if d < (2σ * 1.05)
#         n_cross = -Jij * sin(orientaions[i] - orientaions[j])
#     else
#         n_cross = 0.0
#     end

#     output.torque[i] += n_cross
#     output.torque[j] -= n_cross

#     return output
# end

