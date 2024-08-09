"""
Calculation of order parameter (mean of velocity) in Viscek model. 
"""
function get_order_parameter(all_ori)
    # sample = 1:100:10_000
    # n = length(sample)
    op_list = zeros(100)
    N = length(all_ori)
    for i in 1:100
        op_list[i] = norm(mean(angle2dir.(all_ori[N-i])))
    end
    return mean(op_list)
end

