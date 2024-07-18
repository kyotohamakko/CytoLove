classdef length_ori < handle
  properties
    length
    ori
  end
  methods
    function h = length_ori(length, ori)
      h.length = length;
      h.ori = ori;
    end
  end
end