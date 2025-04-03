classdef node_p < handle
  properties
    point
    next_p
    beg
    ter
  end
  methods
    function h = node_p(point)
      h.point = point;
      h.next_p = -1;
      h.beg = 0.0;
      h.ter = 0.0;
    end
  end
end