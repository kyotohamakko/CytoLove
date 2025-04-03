classdef head_node < handle
  properties
    len_max
    len
    head
  end
  methods
    function h = head_node()
      h.len_max = 2;
      h.len = 0;
      h.head = -1;
    end
    
    function dele(this, target)
        if this.isempty()
            disp("Error in head_node.dele: Try to delete from an empty list! The operation has been retracted.");
            return;
        end
        tem1 = this.head;
        tem2 = this.head.next_p;
        if target == this.head.point
            this.head = tem2;
            this.len = this.len - 1;
            return;
        end
        while isa(tem2,'node')
            if target == tem2.point
                tem1.next_p = tem2.next_p;
                delete(tem2);
                this.len = this.len - 1;
                return;
            else
                tem1 = tem1.next_p;
                tem2 = tem2.next_p;
            end
        end
        if ~isa(tem2,'node')
            disp("Error in head_node.dele: Target not found! The operation has been retracted.");
            disp(tem2);
            return;
        end
    end
    
    function append(this, target)
%         if this.isfull()
%             disp("Error in head_node_p.append: Try to append in full list! The operation has been retracted.");
%             return;
%         end
        
        new_add = node(target);
        new_add.next_p = this.head;
        this.head = new_add;
        this.len = this.len + 1;
    end
    
    function append_norepet(this, target)
        tem = this.head;
        while isa(tem,'node')
            if tem.point == target
                return;
            end
            tem = tem.next_p;
        end
        this.append(target);
    end
    
    function res = isfull(this)
        if this.len == this.len_max
            res=true;
        else
            res=false;
        end
    end
    
    function res = isempty(this)
        if this.len == 0
            res=true;
        else
            res=false;
        end
    end
    
  end
end