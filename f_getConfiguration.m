function [ config_Pa ] = f_getConfiguration( n_levels, nPa )
    config_Pa = 0:(n_levels-1);
    for i = 1 : (nPa-1)
        config_Pa = combvec(config_Pa,0:(n_levels-1));
    end

end

