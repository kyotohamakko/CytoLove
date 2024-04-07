classdef point_stack < handle
  properties
    top
    contant
  end
  methods
    function h = point_stack(contant_datatype_init)
      h.top = 0;
      h.contant = contant_datatype_init;
    end
    function res = isempty(this)
      if this.top == 0
          res = true;
      else
          res = false;
      end
    end
    function push(this, item)
        this.contant(this.top+1) = item;
        this.top = this.top+1;
    end
    function res = pop(this)
        if this.isempty()
            disp("Error in pop in point_stack: Try to pop an empty stack!");
            return;
        end
        res = this.contant(this.top);
        this.contant(this.top) = [];
        this.top = this.top-1;
    end
  end
end