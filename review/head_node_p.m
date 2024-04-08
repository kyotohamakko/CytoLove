classdef head_node_p < handle
  properties
    len_max
    side_len_max
    len
    head
  end
  methods
    function h = head_node_p()
      h.len_max = 2;
      h.side_len_max = 1;
      h.len = 0;
      h.head = -1;
    end
    
    function dele(this, target)
        if this.isempty()
            disp("Error in head_node_p.dele: Try to delete from an empty list! The operation has been retracted.");
            return;
        end
        tem1 = this.head;
        tem2 = this.head.next_p;
        if target == this.head.point
            this.head = tem2;
%             res = this.head.point;
            this.len = this.len - 1;
            return;
        end
        while isa(tem2,'node_p')
            if target == tem2.point
                tem1.next_p = tem2.next_p;
%                 res = tem2;
                delete(tem2);
                this.len = this.len - 1;
                return;
            else
                tem1 = tem1.next_p;
                tem2 = tem2.next_p;
            end
        end
        if ~isa(tem2,'node_p')
            disp("Error in head_node_p.dele: Target not found! The operation has been retracted.");
%             res = -1;
            return;
        end
    end
    
    function append(this, target, beg, ter)
        if this.isfull()
            disp("Error in head_node_p.append: Try to append in full list! The operation has been retracted.");
            return;
        end
        if this.side_full(beg)
            disp("Error in head_node_p.append: Side is full! The operation has been retracted.");
            return;
        end
        new_add = node_p(target);
        new_add.beg = beg;
        new_add.ter = ter;
        new_add.next_p = this.head;
        this.head = new_add;
        this.len = this.len + 1;
    end
    
    function append_norepet(this, target, beg, ter)
        tem = this.head;
        while isa(tem,'node_p')
            if tem.point == target
                return;
            end
            tem = tem.next_p;
        end
        this.append(target, beg, ter);
    end
    
    function res = get(this, idx)
        if idx>this.len
            disp("Error in head_node_p.get: Idx out of range!");
            res = -1;
            return;
        end
        
        if idx == 1
            res = this.head;
            return;
        end
        
        tem = this.head;
        for i=1:1:idx-1
            tem = tem.next_p;
        end
        res = tem;
    end
    
    function res = get_side_conn(this, side, clear_empty_edge)
        tem = this.head;
        
        res = (node_p(point_([-1, -1], 1, 2)));
        if size(res, 2)>this.side_len_max
            disp("Error in head_node_p.get_side_conn: side_len_max changed!");
            return;
        end
        
        i = 1;
        while isa(tem,'node_p')
            if tem.beg == side
                res(i) = tem;
                tem = tem.next_p;
                i=i+1;
            else
                tem = tem.next_p;
            end
        end
        
        if clear_empty_edge && i<=size(res, 2)
            for jj = i:1:size(res, 2)
                res(end) = [];
            end
        end
    end
    
    function res = side_full(this, side)
        if side~=-1.0 && side~=1.0
            disp("Error in head_node_p.side_full: Input should be -1 or 1!");
            return;
        end
        
        res=false;
        tem = this.head;
        count = 0;
        while isa(tem,'node_p')
            if tem.beg == side
                count = count+1;
            end
            tem = tem.next_p;
        end
        
        if count>=this.side_len_max
            res = true;
        end
    end
    
    function res = get_side_len(this, side)
        if side~=-1.0 && side~=1.0
            disp("Error in head_node_p.side_full: Input should be -1 or 1!");
            return;
        end
        
        tem = this.head;
        res = 0;
        while isa(tem,'node_p')
            if tem.beg == side
                res = res+1;
            end
            tem = tem.next_p;
        end
    end
    
    function res = get_empty_side(this)
        res = [-1, 1];
        
        if this.get_side_len(1) ~= 0
            res(2) = [];
        end
        
        if this.get_side_len(-1) ~= 0
            res(1) = [];
        end
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