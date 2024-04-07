classdef total_energy < handle
  properties
    E
    E_array
    E_int_array
    E_ext_array
    p_num_array
  end
  methods
    function h = total_energy(init)
      h.E = init;
      h.E_array(1) = 0.0;
      h.E_int_array(1) = 0.0;
      h.E_ext_array(1) = 0.0;
      h.p_num_array(1) = 0.0;
    end
  end
end