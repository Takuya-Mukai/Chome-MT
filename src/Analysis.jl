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







