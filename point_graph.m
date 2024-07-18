classdef point_graph < handle
  properties
    head
    tail
    ter_tail
    Points
    ter_Points
    pos_map
    point_max_scale
    point_min_scale
  end
  
  methods
    function h = point_graph(map_size, p_scale_max, p_scale_min)
      h.head = 1;
      h.tail = 1;
      h.ter_tail = 1;
      h.point_max_scale = p_scale_max;
      h.point_min_scale = p_scale_min;
      
      point_arr(1) = point_([-1, -1], 1, 2);
      h.Points = point_arr;
      
      ter_point_arr(1) = point_([-1, -1], 1, 2);
      h.ter_Points = ter_point_arr;
      
      for i=1:1:(map_size(1)*map_size(2))
            hnp = head_node();
            p_map(i) = hnp;
      end
      p_map = reshape(p_map, map_size(1), map_size(2));
      h.pos_map = p_map;
    end
    
    function append(this, point)
        this.Points(this.tail) = point;
        this.ter_Points(this.ter_tail) = point;
        this.pos_map(point.pos(2), point.pos(1)).append_norepet(point);
        this.tail = this.tail + 1;
        this.ter_tail = this.ter_tail + 1;
    end
    
    function dele(this, idx)
        if this.isempty()
            disp("Error in graph.dele: Try to delete from an empty list! The operation has been retracted.");
            return;
        end
        if idx>=this.tail
            disp("Error in graph.dele: Index out of range! The operation has been retracted.");
            return;
        end
        target = this.Points(idx);
        this.pos_map(target.pos(2), target.pos(1)).dele(target);
        
        res = find(this.ter_Points==target);
        if size(res, 2)==0 || size(res, 2)>1
            disp("Error in graph.dele: Unknow Error!");
            return;
        end
        
        delete(target);
        
        this.Points(idx) = [];
        this.tail = this.tail - 1;
        
        this.ter_Points(res) = [];
        this.ter_tail = this.ter_tail - 1;
    end
    
    function dele_(this, target, ter_id)
        if this.isempty()
            disp("Error in graph.dele: Try to delete from an empty list! The operation has been retracted.");
            return;
        end
        if ~ismember(target, this.Points)
            disp("Error in graph.dele: Try to delete unexist target! The operation has been retracted.");
            return;
        end
        
        this.pos_map(target.pos(2), target.pos(1)).dele(target);
        
        res = find(this.Points==target);
        if size(res, 2)==0 || size(res, 2)>1
            disp("Error in graph.dele: Unknow Error!");
            return;
        end
        
        delete(target);
        
        this.Points(res) = [];
        this.tail = this.tail - 1;
        
        if ter_id>0
            this.ter_Points(ter_id) = [];
            this.ter_tail = this.ter_tail - 1;
        end
    end
    
    function dele_conn(this, p1, p2)
        if ismember(p1, this.Points) && ismember(p2, this.Points)
            p1.head_node_p.dele(p2);
            p2.head_node_p.dele(p1);
        else
            disp("Error in point_graph.dele_conn: Try to dele conn that does't exist!");
            return;
        end
        
        this.ter_adjust(p1);
        this.ter_adjust(p2);
    end
    
    function res = add_conn(this, p1, conn, p1_id)
        p2 = conn.point;
        p1_beg = conn.beg;
        p1_ter = conn.ter;
        
        res = false;
        
        if p2.pos(1)<0 || p2.pos(2)<0
            return;
        end
        
        if ismember(p1, this.Points) && ismember(p2, this.Points)
            p1.head_node_p.append_norepet(p2, p1_beg, p1_ter);
            p2.head_node_p.append_norepet(p1, p1_ter, p1_beg);
        else
            disp("Error in point_graph.dele_conn: Try to dele conn that does't exist!");
            return;
        end
        
        if p1_id>0
            if this.loop_check(p1_id)
                this.dele_conn(p1, conn.point);
            else
                res = true;
            end
        end
        
        this.ter_adjust(p1);
        this.ter_adjust(p2);
    end
    
    function ter_adjust(this, p)
        res = find(this.ter_Points==p);
        
        if p.head_node_p.len<=1 && size(res, 2)<=0
            this.ter_Points(this.ter_tail) = p;
            this.ter_tail = this.ter_tail + 1;
        elseif p.head_node_p.len>1 && size(res, 2)>0
            this.ter_Points(res) = [];
            this.ter_tail = this.ter_tail - 1;
        end
    end
    
    function p_move_rotate_scale(this, p, pos_new, ori_new)
        p_in_newpos = this.pos_map(pos_new(2), pos_new(1));
        if p_in_newpos.len>0
            return;
        end
        
        old_pos = p.pos;
%         old_orient = p.orient;
        
        p.pos = pos_new;
        p.orient = ori_new;
        
        this.pos_map(old_pos(2), old_pos(1)).dele(p);
        this.pos_map(pos_new(2), pos_new(1)).append_norepet(p);
        
%         if this.coll_check(p, old_pos)
%             this.pos_map(old_pos(2), old_pos(1)).dele(p);
%             this.pos_map(pos_new(2), pos_new(1)).append_norepet(p, -1, -1);
%         else
%             p.pos = old_pos;
%             p.orient = old_orient;
%         end
    end
    
    function res = loop_check(this, p_idx)
        res = false;

        visit = false*(ones(size(this.Points)));
        pull = false*(ones(size(this.Points)));
        para = -1*(ones(size(this.Points)));

        ps = point_stack(1);

        ps.push(p_idx);
        pull(p_idx)=true;

        while ~ps.isempty()
            item_curr = ps.pop();
            p_curr = this.Points(item_curr);
            visit(item_curr) = true;
            nl = p_curr.head_node_p;
            if ~nl.isempty()
                tem = nl.head;
                while isa(tem,'node_p')
                    idx = find(this.Points==tem.point);
                    para(idx)=item_curr;

                    if ~visit(idx) && ~pull(idx)
                        ps.push(idx);
                        pull(idx)=true;
                    end
                    if visit(idx) && para(item_curr) ~= idx
                        res = true;
                        return;
                    end
                    tem = tem.next_p;
                end
            end
        end
    end
    
    function res=coll_check(this, p, skip, cw1, cw2)
        res=0.0;
        shift = 2*this.point_max_scale;
        s = size(this.pos_map);
        p_pos = p.pos;
        for i=(p_pos(2)-shift):1:(p_pos(2)+shift)
            for j=(p_pos(1)-shift):1:(p_pos(1)+shift)
                if i>0 && j>0 && i<=s(1) && j<=s(2)
                    
                    tem_hnp = this.pos_map(i, j).head;
                    while isa(tem_hnp,'node')
                        tem_p = tem_hnp.point;
                        if tem_p.pos(1) > 0 && tem_p.pos(2) > 0 && tem_p~=skip
                            if this.intersect(tem_p, p)
                                E1 = (tem_p.scale+p.scale)/2;
                                E2 = cos(pi*(tem_p.orient-p.orient)/180)^(2*cw2);
                                E = E1 * (1+cw1*E2);
                                res = res+E;
                            end
                        end
                        tem_hnp = tem_hnp.next_p;
                    end
                    
                end
            end
        end
    end

    function res = intersect(this, p1, p2)
        this.length;
        cir1 = [p1.pos(1), p1.pos(2), p1.scale, p1.scale/4.0, p1.orient+90.0];
        cir2 = [p2.pos(1), p2.pos(2), p2.scale, p2.scale/4.0, p2.orient+90.0];

        fcir = @(cir)this.calculateEllipse(cir(1),cir(2),cir(3),cir(4),cir(5), 100);
        
        [x1,y1] = fcir(cir1);
        [x2,y2] = fcir(cir2);
        
        if p1.scale > p2.scale
            idx = this.incir(cir1,x2,y2);
        else
            idx = this.incir(cir2,x1,y1);
        end
        
        if sum(idx) > 0
            res = true;
        else
            res = false;
        end
    end
    
    function F = incir(this, cir,x,y)
        this.length;
        x = x-cir(1);
        y = y -cir(2);
        cir(5) = cir(5)*pi/180;
        T =[cos(cir(5)), sin(cir(5));-sin(cir(5)),cos(cir(5))];
        xy = T*[x';y'];
        xy(1,:) = xy(1,:)+cir(1) ;
        xy(2,:) = xy(2,:)+cir(2) ;
        F =((xy(1,:)'-cir(1)).^2/cir(3)^2 + (xy(2,:)'-cir(2)).^2/cir(4)^2)<=1;
    end
    
    function [X,Y] = calculateEllipse(this, x, y, a, b, angle, steps)
        this.length;

        beta = angle * (pi / 180);
        sinbeta = sin(beta);
        cosbeta = cos(beta);

        alpha = linspace(0, 360, steps)' .* (pi / 180);
        sinalpha = sin(alpha);
        cosalpha = cos(alpha);

        X = x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
        Y = y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);

        if nargout==1, X = [X Y]; end
    end
    
    function res = isempty(this)
        if this.tail == this.head
            res=true;
        else
            res=false;
        end
    end
    
    function res = is_ter_empty(this)
        if this.ter_tail == this.head
            res=true;
        else
            res=false;
        end
    end
    
    function res = length(this)
        res = this.tail - this.head;
    end
    
    function plot_all_points(this)
        for i=1:1:size(this.Points, 2)
            p = this.Points(i);
            pos = p.pos;
            orient = p.orient;
            scale = p.scale;
            birth_type = p.birth_type;
            
            [ex,ey] = this.calculateEllipse(pos(1), pos(2), scale, scale/4.0, orient+90.0, 100);
            
            if birth_type == 1
                plot(pos(1), pos(2),'.r','MarkerSize',5);
                plot(ex,ey, 'LineWidth',2, 'Color',[1 0 0]);
            elseif birth_type == 2
                plot(pos(1), pos(2),'.g','MarkerSize',5);
                plot(ex,ey, 'LineWidth',2, 'Color',[0 1 0]);
            elseif birth_type == 3
                plot(pos(1), pos(2),'.k', 'MarkerSize',5);
                plot(ex,ey, 'LineWidth',2, 'Color',[0.4940 0.1840 0.5560]);
            end
        end
    end
    
    function plot_ter_points(this)
        for i=1:1:size(this.ter_Points, 2)
            p = this.ter_Points(i);
            pos = p.pos;
            orient = p.orient;
            scale = p.scale;
            birth_type = p.birth_type;
            
            [ex,ey] = this.calculateEllipse(pos(1), pos(2), scale, scale/4.0, orient+90.0, 100);
            if birth_type == 1
                plot(pos(1), pos(2),'.r','MarkerSize',5);
                plot(ex,ey, 'LineWidth',2, 'Color',[1 0 0]);
            elseif birth_type == 2
                plot(pos(1), pos(2),'.g','MarkerSize',5);
                plot(ex,ey, 'LineWidth',2, 'Color',[0 1 0]);
            end
        end
    end
    
    function plot_connect_pos(this)
        for i=1:1:size(this.Points, 2)
            nl = this.Points(i).head_node_p;
            if nl.isempty()
                continue;
            end
            pos_father = this.Points(i).pos;
            tem = nl.head;
            while isa(tem,'node_p')
                pos_child = tem.point.pos;
                line([pos_father(1) pos_child(1)],[pos_father(2) pos_child(2)],'color','b','linewidth',3);
                tem = tem.next_p;
            end
        end
    end
    
    function plot_connect(this)
        for i=1:1:size(this.Points, 2)
            p1 = this.Points(i);
            nl = p1.head_node_p;
            if nl.isempty()
                continue;
            end
            
            pos1 = p1.pos;
            orient1 = p1.orient;
            scale1 = p1.scale;
            n1 = (scale1/4.0).*[cos((orient1/180)*pi), sin((orient1/180)*pi)];
            
            tem = nl.head;
            while isa(tem,'node_p')
                p2 = tem.point;
                pos2 = p2.pos;
                orient2 = p2.orient;
                scale2 = p2.scale;
                
                beg = tem.beg;
                ter = tem.ter;
                
                n2 = (scale2/4.0).*[cos((orient2/180)*pi), sin((orient2/180)*pi)];
                pos1_ = pos1 + beg*n1;
                pos2 = pos2 + ter*n2;
                
%                 plot(pos1(1), pos1(2),'.r','MarkerSize',5);
%                 plot(pos2(1), pos2(2),'.r','MarkerSize',5);
                
                line([pos1_(1) pos2(1)],[pos1_(2) pos2(2)],'color','b','linewidth',3);
                tem = tem.next_p;
            end
        end
    end
    
  end
end