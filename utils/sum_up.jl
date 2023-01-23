using DataFrames;


function sum_up_lshape(res_l_shape)
    
    x = pop!(res_l_shape, "first_level_var", nothing)
    d = DataFrame(res_l_shape)

    return x, d
end;