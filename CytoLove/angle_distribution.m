clear
clc

actin_type = "lamilipodia";
load(strcat("save/", actin_type, "_graph.mat"),"pg");

pg.Points = pg.Points(randperm(size(pg.Points,2)));

figure;imshow(ones(size(pg.pos_map)));title("plot in different color");
hold on;

colors = ['y','m','c','r','g','b','w','k'];
ploted_idx = zeros(size(pg.Points));


%%%%%%%%%%%%%%%%%%  Mark particles in shot actin  %%%%%%%%%%%%%%%%%%
idx = 1;
len_thre = 10;
for i=1:1:size(pg.Points, 2)
    if ploted_idx(i)==0
        ploted_idx =mark_noplot(pg, i, ploted_idx, len_thre);
        idx = idx+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


stat = length_ori(-1, -1);

idx = 1;
for i=1:1:size(pg.Points, 2)
    if ploted_idx(i)==0
        col_idx = mod(idx, size(colors, 2))+1;
        [ploted_idx, stat] =graph_analysis(pg, i, ploted_idx, stat, false, colors(col_idx));
        idx = idx+1;
    end
end

[~, ind] = sort([stat.ori]);
stat_sort = stat(ind);


%%%%%%%%%%%%%%%%%%  Angle Distribution  %%%%%%%%%%%%%%%%%%
bin_len = 7.5;
bin_star = bin_len;
total_len = 0.0;
part_len = 0.0;

stat_res = zeros(1, ceil(180/bin_len));
xlable = zeros(size(stat_res));%[15,30,45,60,75,90,105,120,135,150,165,180];
for lable_id=1:1:size(xlable, 2)
    if lable_id == 1
        xlable(lable_id) = 0.0-90.0;
    else
        xlable(lable_id) = xlable(lable_id-1)+bin_len;
    end
end

idx_stat_res = 1;

for idx_stat=1:1:size(stat_sort, 2)
    item = stat_sort(idx_stat);
    ori_ = item.ori;
    length_ = item.length;
    
    total_len = total_len+length_;
    
    if ori_>=180
        bin_id = size(stat_res, 2);
    elseif mod(ori_, bin_len) == 0
        bin_id = ceil(ori_/bin_len)+1;
    else
        bin_id = ceil(ori_/bin_len);
    end
    stat_res(bin_id) = stat_res(bin_id)+length_;
end

stat_res = stat_res/total_len;

figure;bar(xlable, stat_res, 'histc');
ylim([0 0.1]);
yticks([0 0.025 0.05 0.075 0.1])
xticks([-90 -45 0 45 90])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function res1 = mark_noplot(pg, p_idx, ploted_idx, len_thre)
    visit = false*(ones(size(pg.Points)));
    pull = false*(ones(size(pg.Points)));

    ps = point_stack(1);

    ps.push(p_idx);
    pull(p_idx)=true;
    
    
    sgs_idx = 1;
    while ~ps.isempty()
        item_curr = ps.pop();
        p_curr = pg.Points(item_curr);
        
        sub_graph_stat(sgs_idx) = item_curr;
        sgs_idx = sgs_idx+1;
        
        visit(item_curr) = true;
        nl = p_curr.head_node_p;
        
        if ~nl.isempty()
            tem = nl.head;
            while isa(tem,'node_p')
                idx = find(pg.Points==tem.point);

                if ~visit(idx) && ~pull(idx)
                    ps.push(idx);
                    pull(idx)=true;
                end
                tem = tem.next_p;
            end
        end
    end
    
    if size(sub_graph_stat, 2)<=len_thre
        ploted_idx(sub_graph_stat) = 1;
    end
    
    res1 = ploted_idx;
end


function [res1, res2] = graph_analysis(pg, p_idx, ploted_idx, stat, plot_point, col)
    col=rand(1,3);
    visit = false*(ones(size(pg.Points)));
    pull = false*(ones(size(pg.Points)));

    ps = point_stack(1);

    ps.push(p_idx);
    pull(p_idx)=true;
    
    ploted_idx(p_idx) = 1;

    while ~ps.isempty()
        item_curr = ps.pop();
        p_curr = pg.Points(item_curr);
        
        ploted_idx(item_curr) = 1;
        
        if plot_point
            pos = p_curr.pos;
            orient = p_curr.orient;
            scale = p_curr.scale;

            [ex,ey] = calculateEllipse(pos(1), pos(2), scale, scale/4.0, orient+90.0, 100);
            plot(pos(1), pos(2),'.r','MarkerSize',5);
            plot(ex,ey, 'LineWidth',2, 'Color',[1 0 0]);
        end
        
        
        visit(item_curr) = true;
        nl = p_curr.head_node_p;
        
        if ~nl.isempty()
            pos_father = p_curr.pos;
            
            tem = nl.head;
            while isa(tem,'node_p')
                idx = find(pg.Points==tem.point);
                if ploted_idx(idx) == 1
                    tem = tem.next_p;
                    continue;
                else
                    ploted_idx(idx) = 1;
                end
                
                pos_child = tem.point.pos;
                line([pos_father(1) pos_child(1)],[pos_father(2) pos_child(2)],'color',col,'linewidth',3);
                
                lengh = sqrt( (pos_father(1)-pos_child(1))^2 + (pos_father(2)-pos_child(2))^2 );
                
                ori = (1/2)*(p_curr.orient+tem.point.orient);
                
                if ori<0
                    ori = ori+180.0;
                end
                
                if stat(end).length ~= -1
                    stat(end+1) = length_ori(lengh, ori);
                else
                    stat(end) = length_ori(lengh, ori);
                end

                if ~visit(idx) && ~pull(idx)
                    ps.push(idx);
                    pull(idx)=true;
                end
                tem = tem.next_p;
            end
        end
    end
    
    res1 = ploted_idx;
    res2 = stat;
end

function [X,Y] = calculateEllipse(x, y, a, b, angle, steps)

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