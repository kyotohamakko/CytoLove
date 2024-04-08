classdef data_load < handle
  properties
    AODF_F
    theta
  end
  methods
    function h = data_load(d)
      h.AODF_F = d;
      
      L = size(d, 1);
      id = 1:L;
      theta_ = angle(exp(2*pi*1i*id/L))*180/pi;
      theta_ = theta_(theta_>0);
      h.theta = theta_;
    end
  end
end