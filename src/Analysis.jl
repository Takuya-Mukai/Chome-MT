"""
Calculation of order parameter (mean of velocity) in Viscek model. 
"""
function get_order_parameter(all_ori)
    op_list = zeros(10)
    N = length(all_ori)
    for i in 1:10
        op_list[i] = norm(mean(angle2dir.(all_ori[N-i])))
    end
    return mean(op_list)
end


function get_XYorder_parameter(all_ori)
    x_list = zeros(10)
    y_list = zeros(10)
    N = length(all_ori)
    for i in 1:10
        x_list[i], y_list[i] = mean(angle2dir.(all_ori[N-i]))
    end
    return mean(abs.(x_list)), mean(abs.(y_list))
end

function is_inside_region(pos, xmin, xmax, ymin, ymax)
    x, y = pos
    return (xmin <= x <= xmax) && (ymin <= y <= ymax)
end


function get_order_parameter_inside(all_pos, all_ori, xmin, xmax, ymin, ymax)
    op_list = zeros(10)
    N = length(all_ori)
    for i in 1:10
        p_inside = all_ori[N-i][is_inside_region.(all_pos[N-i], xmin, xmax, ymin, ymax)]
        op_list[i] = norm(mean(angle2dir.(p_inside)))
    end
    return mean(op_list)
end

function get_XYorder_parameter_inside(all_pos, all_ori, xmin, xmax, ymin, ymax)
    x_list = zeros(10)
    y_list = zeros(10)
    N = length(all_ori)
    # Np = length(all_ori[end])
    for i in 1:10
        p_inside = all_ori[N-i][is_inside_region.(all_pos[N-i], xmin, xmax, ymin, ymax)]
        x_list[i], y_list[i] = mean(angle2dir.(p_inside))
    end
    return mean(abs.(x_list)), mean(abs.(y_list))
end




