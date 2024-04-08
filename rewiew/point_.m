classdef point_ < handle
  properties
    pos
    orient
    scale
    birth_type
    head_node_p
  end
  methods
    function h = point_(pos, s_min, s_max)
      h.pos = pos;
      h.orient = unifrnd(0, 180);
      
      if s_min==s_max
        h.scale = s_max;
      else
          h.scale = unifrnd(s_min, s_max);
      end
      
      h.birth_type = -1;
      h.head_node_p = head_node_p();
    end
  end
end