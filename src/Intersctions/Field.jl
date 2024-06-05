"""
    Updata field 
"""

function UpdateField(pos, field, field0, paras, dt)
    σ=paras.sigma/6
    dx = paras.dx
    A = paras.A 
    λ = .0
    L = paras.L
    xref = SA[L/2 , L/2]
    sides = SA[L , L]

    #* decay
    # field .-= field0 .* λ*dt
    field .= field0 .* exp(-λ*dt)

    #* add source
    for i in eachindex(pos)
        r0 = pos[i]
        lo = r0 .- 4σ
        up = r0 .+ 4σ
        #* mapping physical position to grid index
        lo_idx, up_idx = pos2idx(lo, dx), pos2idx(up, dx)
        @inbounds Threads.@threads for yj in lo_idx[2]:up_idx[2]
            for xi in lo_idx[1]:up_idx[1]
                r = idx2pos(SA[xi, yj], dx)
                # wrapped_pos = wrap_relative_to(r, xref, sides)
                wrapped_pos = wrap_position.(r, sides)
                xii, yjj = pos2idx(wrapped_pos, dx)
                field[xii, yjj] += distribute(r, r0, σ, A) * dt
            end
        end
    end

    field[field .> 1.5] .= 1.5
    # updataGrid!(field, field0, paras, dt)
    # PeriodicBoundary!(field0, field, paras, dt)
    
    return field, field0
end



function distribute(pos, pos0, σ, A)
    dr = norm(pos - pos0)
    return A* exp(-(dr)^2 / (4 * σ )) / (4π * σ )
end


# exp_decay(λ ,dt) = exp(-λ * dt)
# linear_decay(λ, dt) = 



#= 
* update grid in diffusion equation with 2-order method
=#
function updataGrid!(u, du, paras, dt)
    dx = paras.dx
    dy = dx
    nx, ny = size(u)
    D = 1/50
    factor = dt * D
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    Threads.@threads for j in 2:ny-1
        for i in 2:nx-1
            @inbounds du[i, j] = u[i, j] + factor * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2
                                                     +
                                                     (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2)
        end
    end
end


function PeriodicBoundary!(u, du, paras, dt)
    dx = paras.dx
    dy = dx
    nx, ny = size(u)
    D = 0.05
    factor = dt * D
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    #* updata the top and bottom along x
    for i in 2:nx-1
        du[i, 1] = u[i, 1] + dt * D * ((u[i+1, 1] - 2 * u[i, 1] + u[i-1, 1]) / dx^2
                                    +
                                    (u[i, 2] - 2 * u[i, 1] + u[i, ny]) / dy^2)

        du[i, ny] = u[i, ny] + dt * D * ((u[i+1, ny] - 2 * u[i, ny] + u[i-1, ny]) / dx^2
                                        +
                                        (u[i, 1] - 2 * u[i, ny] + u[i, ny-1]) / dy^2)
    end

    #* updata the left and right along y
    for j in 2:ny-1
            du[1, j] = u[1, j] + dt * D * ((u[2, j] - 2 * u[1, j] + u[nx, j]) / dx^2 
                                           +
                                           (u[1, j+1] - 2 * u[1, j] + u[1, j-1]) / dy^2)

            du[nx, j] = u[nx, j] + dt * D * ((u[1, j] - 2 * u[nx, j] + u[nx-1, j]) / dx^2 
                                            +
                                            (u[nx, j+1] - 2 * u[nx, j] + u[nx, j-1]) / dy^2)                              
    end

    #* updata 4 corner
    du[1,1] = u[1,1] + dt * D *((u[2,1] - 2*u[1,1] + u[nx,1]) / dx^2 + (u[1,2] - 2*u[1,1] + u[1,ny]) / dy^2)
    du[nx,ny] = u[nx,ny] + dt * D *((u[1,ny] - 2*u[nx,ny] + u[nx-1,ny]) / dx^2 + (u[nx,1] - 2*u[nx,ny] + u[nx,ny-1]) / dy^2)
    du[1,ny] = u[1,ny] + dt * D *((u[2,ny] - 2*u[1,ny] + u[nx,ny]) / dx^2 + (u[1,1] - 2*u[1,ny] + u[1,ny-1]) / dy^2)
    du[nx,1] = u[nx,1] + dt * D *((u[1,1] - 2*u[nx,1] + u[nx-1,1]) / dx^2 + (u[nx,2] - 2*u[nx,1] + u[nx,ny]) / dy^2)
end
