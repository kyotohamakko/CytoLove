classdef node < handle
  properties
    point
    next_p
  end
  methods
    function h = node(point)
      h.point = point;
      h.next_p = -1;
    end
  end
end