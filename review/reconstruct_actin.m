clear
clc

actin_type = "cortex";
src = im2double(imread(strcat("data/", actin_type, "_frame.bmp")));

%%%%%%%%%%% Parameter %%%%%%%%%%%
if actin_type == "cortex"
    wid = 3.0;
elseif actin_type == "lamilipodia"
    wid = 2.0;
end

object_remove = 20;

lam = 0.1/4;

p_scale_max = 2.0;
p_scale_min = 2.0;

k = 2.0;
c = -4.0;
c_m = -15.0;
k_track = 0.2;
T = 1.0;

elmax = floor(2*p_scale_max);
bta1 = 2.0;
bta2 = 1.0;
T_conn = 1.0;
b = 2.0;

cw1 = 100.0;
cw2 = 6.0;

TL = 0.05;
TH = 0.15;

typ = 1;
con_type = -1;

iter=5000000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ei, orient, AODF_F, ei2, orient2] = sdeconv(src,'reg',0.5,'wid',wid,'lam', lam, 'aniso',1.0,'numits',50);
ei(ei<=0)=0.0;
ei = imadjust(ei,stretchlim(ei),[]);

[src_tem, ei, orient, AODF_F] = res_trim(src, ei, orient, AODF_F, 16);
orient_test = angle(orient)*180/pi;

mask = pos_mask(ei, orient, object_remove, TL, TH);
[ei, AODF_F] = data_masked(ei, AODF_F, mask);

[p_birth, ei_nor] = ei_normalization(ei, 2.0*sqrt(wid), mask);

pos = 1:1:size(p_birth, 1);

proposal = ["birth", "death", "connect", "move_rotate", "track_birth", "track_death"];%
pg = point_graph(size(ei), p_scale_max, p_scale_min);

dl = data_load(AODF_F);

E_total = total_energy(0.0);

disp("Start...");

tic
for i=1:1:iter
    if mod(i, 10000)==0
        disp(i);
    end
    
    pick = proposal(randsample(1:1:size(proposal, 2),1));
    
    switch pick
       case 'birth'
          birth(pg, ei, orient_test, pos, p_scale_min, p_scale_max, p_birth, cw1, cw2, k, c, T, E_total);
       case 'death'
          death(pg, ei, p_birth, cw1, cw2, k, c, T, E_total);
       case 'connect'
          connect(pg, elmax, bta1, bta2, T_conn, b, T, E_total, con_type);
       case 'move_rotate'
          move_rotate_rescale(pg, dl, ei, c_m, T, bta1, bta2, cw1, cw2, E_total);
       case 'track_birth'
           track_birth(pg, dl, ei, orient_test, p_scale_min, p_scale_max, bta1, bta2, cw1, cw2, k_track, c_m, T, typ, E_total);
       case 'track_death'
           track_death(pg, dl, bta1, bta2, cw1, cw2, k_track, c_m, T, E_total);
       otherwise
          disp('Unknow error!')
    end
    
    if mod(i, 20000)==0
        T = T*0.9995;
    end

    if mod(i, 1000000) == 0
            typ = 2;
    end

    if mod(i, 4000000) == 0
        con_type = 2;
        proposal = ["connect", "track_birth", "track_death"];
    end
end
toc

save(strcat("save/", actin_type, "_graph.mat"));

figure;imshow(src_tem);
hold on;

pg.plot_all_points();
pg.plot_connect_pos();

figure;
yyaxis left
plot(E_total.E_array, 'LineWidth',2, 'color', [1 0 0]);
hold on;
plot(E_total.E_int_array, 'LineWidth',2, 'color', [0 1 0]);
hold on;
plot(E_total.E_ext_array, 'LineWidth',2, 'color', [0 0 1]);
hold on;
ylabel('Energy', 'FontSize', 20);
yyaxis right
plot(E_total.p_num_array, 'LineWidth',2, 'color', [0 0 0]);
legend('E\{X\}(M)','E_{int}(M)','E_{ext}(M,X)','number of particles', 'location', 'east');
xlabel('Iteration', 'FontSize', 20);
ylabel('Number of particles', 'FontSize', 20);

function [res1, res2] = data_masked(ei, AODF_F, mask)
    res1 = ei.*mask;
    res2 = ones(size(AODF_F));
    for i = 1:1:size(AODF_F, 1)
        res2(i,:,:) = squeeze(AODF_F(i,:,:)).*mask;
    end
end

function [res, ei_nor] = ei_normalization(ei, sig, mask)
    gauss = single(circshift(fspecial('gaussian',size(ei),sig),[1 1]));
    ftgauss = fft2(fftshift(gauss));

    mean = (real(ifft2(fft2(ei).*ftgauss)));
    var = (real(ifft2(fft2(ei.^2).*ftgauss))) - mean.^2;
    var(var<=0)=0.0;
    
    ei_nor = (ei-mean)./(sqrt(var)+0.001);
    
    ei_nor(ei_nor<=0.0) = 0.0;
    ei_nor = ei_nor.*mask;
    ei_nor_T = (ei_nor/sum(sum(ei_nor)))';
    res = ei_nor_T(:);
end

function [src_tem, ei_, ori, AODF] = res_trim(src, ei, orient, AODF_F, sub)
    src_tem = src(sub:end-sub, sub:end-sub);
    ei_ = ei(sub:end-sub, sub:end-sub);
    ori = orient(sub:end-sub, sub:end-sub);
    AODF = AODF_F(:, sub:end-sub, sub:end-sub);
end

function [src_tem, ei_, ori, AODF] = res_trim2(src, ei, orient, AODF_F, ver, hor)
    src_tem = src(ver(1):ver(2), hor(1):hor(2));
    ei_ = ei(ver(1):ver(2), hor(1):hor(2));
    ori = orient(ver(1):ver(2), hor(1):hor(2));
    AODF = AODF_F(:, ver(1):ver(2), hor(1):hor(2));
end

function res = pos_mask(ei, orient, object_remove, TL, TH)
    temp = non_maxsupress(ei, 1j*orient, TL, TH);
    temp(temp>0) = true;
    temp(temp<=0) = false;
    if object_remove>0
        temp = bwareaopen(temp, object_remove);
    end
    B=[0 1 0
       1 1 1
       0 1 0];
    temp=imdilate(temp,B);
    temp = bwmorph(temp,'thin',Inf);
    res = temp;
end

function birth(pg, ei, orient_test, pos, p_scale_min, p_scale_max, p_birth, cw1, cw2, k, c, T, E_total)
    y = randsample(pos,1,true,p_birth);
    [col_idx, row_idx] = pos_convert_stack(y, size(ei, 2));
    p = point_([col_idx, row_idx], p_scale_min, p_scale_max);
    p.orient = orient_test(row_idx, col_idx);
    p.birth_type = 1;
    
    E_coll = pg.coll_check(p, -1, cw1, cw2);
    
    [R_bir, dE] = birth_R(ei(row_idx, col_idx), E_coll, p_birth(y), 1.0/(pg.length()+1), k, c, T);
    if ~reject(R_bir)
        pg.append(p);
        E_total.E = E_total.E + dE;
        E_total.E_array(end+1) = E_total.E;
        E_total.E_int_array(end+1) = E_total.E_int_array(end);
        E_total.E_ext_array(end+1) = E_total.E_ext_array(end) + dE;
        E_total.p_num_array(end+1) = E_total.p_num_array(end) + 1;
    end
end

function [R, dE] = birth_R(A_x, E_coll, P_t, P_t_, k, c, T)
    dE = E_coll + (A_x-k)*c;
    R = exp(-dE/T) * (P_t_/P_t);
end

function death(pg, ei, p_birth, cw1, cw2, k, c, T, E_total)
    if pg.isempty()
        return;
    end
    
    pidx = randsample(1:1:pg.length(),1);
    p = pg.Points(pidx);
    y_death = pos_convert_flat(p.pos, size(ei, 2));
    
    E_coll = pg.coll_check(p, p, cw1, cw2);
    
    [R_dea, dE] = death_R(ei(p.pos(2), p.pos(1)), E_coll, 1.0/pg.length(), p_birth(y_death), k, c, T);
    if ~reject(R_dea) && p.head_node_p.isempty()
        pg.dele(pidx);
        E_total.E = E_total.E + dE;
        E_total.E_array(end+1) = E_total.E;
        E_total.E_int_array(end+1) = E_total.E_int_array(end);
        E_total.E_ext_array(end+1) = E_total.E_ext_array(end) + dE;
        E_total.p_num_array(end+1) = E_total.p_num_array(end) - 1;
    end
end

function [R, dE] = death_R(A_x, E_coll, P_t, P_t_, k, c, T)
    dE = -E_coll - (A_x-k)*c;
    R = exp(-dE/T) * (P_t_/P_t);
end

function connect(pg, elmax, bta1, bta2, T_conn, b, T, E_total, con_type)
    if pg.isempty()
        return;
    end
    
    if con_type <=0
        con_type = randsample([1,2],1,true,[1.0/2, 1.0/2]);
    end
    
    if con_type==1
        arr_len = size(pg.Points, 2);
        p_id = randsample(1:1:arr_len, 1);
        p1 = pg.Points(p_id);
        p_pos = p1.pos;

        side = randsample([-1.0, 1.0],1,true,[1.0/2, 1.0/2]);
    elseif con_type==2
        pidx = randsample(1:1:size(pg.ter_Points, 2),1);
        p1 = pg.ter_Points(pidx);
        p_pos = p1.pos;
        p_id = find(pg.Points==p1);
        
        if size(p_id, 2)==0 || size(p_id, 2)>1
            disp("Error in connect: ox000!");
            return;
        end

        empty_sides = p1.head_node_p.get_empty_side();

        if size(empty_sides, 2) == 0
            disp("Error in connect: ox001!");
            return;
        end

        if size(empty_sides, 2) == 2
            side = randsample(empty_sides,1,true,[1.0/2, 1.0/2]);
        else
            side = empty_sides;
        end
    else
        disp("Error in connect: Unknow connection propose.");
        return;
    end
    
    side_conn = p1.head_node_p.get_side_conn(side, false);
    e0 = side_conn(randsample(1:1:size(side_conn, 2), 1));
    e0 = e_copy(e0);
    p1_adj = e0.point;
    if p1_adj.pos(1)<0 || p1_adj.pos(2)<0
        if p1.head_node_p.isfull()
            return;
        end
    else
        pg.dele_conn(p1, p1_adj);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_p1 = side*(p1.scale/4.0).*[cos((p1.orient/180)*pi), sin((p1.orient/180)*pi)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = size(pg.pos_map);
    
    idx_tem = 1;
    e0_idx = -1;
    
    for i=(p_pos(2)-elmax):1:(p_pos(2)+elmax)
        for j=(p_pos(1)-elmax):1:(p_pos(1)+elmax)
            if i>0 && j>0 && i<=s(1) && j<=s(2) && (i~=p_pos(2) || j~=p_pos(1))
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if dot(n_p1, [j, i]-p1.pos)<0
                    continue;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                tem_hnp = pg.pos_map(i, j).head;
                while isa(tem_hnp,'node')
                    tem_p = tem_hnp.point;

                    if tem_p.pos(1) > 0 && tem_p.pos(2) > 0 && ~tem_p.head_node_p.isfull() %&& ~pg.intersect(tem_p, p1)
                        
                        side2 = closest_side(p1, side, tem_p);
                        if ~tem_p.head_node_p.side_full(side2)
                            e_can = node_p(tem_p);
                            e_can.beg = side;
                            e_can.ter = side2;
                            e_candidate(idx_tem) = e_can;

                            C_edge = C_Edge(p1, side, tem_p, side2, bta1, bta2);
                            p_e(idx_tem) = exp(-(C_edge-b)/T_conn);

                            if sum(tem_p.pos - p1_adj.pos) == 0 && e0.beg == side && e0.ter == side2
                                e0_idx = idx_tem;
                            end

                            idx_tem = idx_tem+1;
                        end
                        
                    end
                    
                    tem_hnp = tem_hnp.next_p;
                
                end

            end
        end
    end

    if ~exist('e_candidate','var') || ~exist('p_e','var')
        return;
    end
    
    p_e(idx_tem) = min(p_e);
    e_candidate(idx_tem) = node_p(point_([-1, -1], 1, 2));

    p_e = p_e./sum(p_e);

%     e_new_id = randsample(1:1:size(e_candidate, 2),1,true,p_e);
    [p_max, e_new_id] = max(p_e);
    e_new = e_candidate(e_new_id);
    
    if e_new.point.pos(1)<0 || e_new.point.pos(2)<0
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if e0_idx<0 && (e0.point.pos(1)>0 || e0.point.pos(2)>0)
        return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [R, E0, dE] = conn_R(p1, e0, e0_idx, e_new, e_new_id, p_e, 2.0, bta2, T);
    
    if ~reject(R) && sum(e_new.point.pos - p1_adj.pos) ~= 0
        succ = pg.add_conn(p1, e_new, p_id);
        if ~succ && (p1_adj.pos(1)>0 || p1_adj.pos(2)>0)
            pg.add_conn(p1, e0, -1);
            E_total.E = E_total.E - E0;
            E_total.E = E_total.E + dE;
            E_total.E_array(end+1) = E_total.E;
            
            E_total.E_int_array(end+1) = E_total.E_int_array(end) - E0 + dE;
            E_total.E_ext_array(end+1) = E_total.E_ext_array(end);
            E_total.p_num_array(end+1) = E_total.p_num_array(end);
        end
    elseif (p1_adj.pos(1)>0 || p1_adj.pos(2)>0)
        pg.add_conn(p1, e0, -1);
    end
end

function res = closest_side(p1, side, tem_p)
    p1_n = (p1.scale/4.0).*[cos((p1.orient/180)*pi), sin((p1.orient/180)*pi)];
    tem_p_n = (tem_p.scale/4.0).*[cos((tem_p.orient/180)*pi), sin((tem_p.orient/180)*pi)];
    
    p1_side = p1.pos + side*p1_n;
    
    tem_p_side_plus = tem_p.pos + tem_p_n;
    tem_p_side_minus = tem_p.pos - tem_p_n;
    
    d_plus = sum((p1_side-tem_p_side_plus).^2);
    d_minus = sum((p1_side-tem_p_side_minus).^2);
    
    if d_plus>d_minus
        res = -1.0;
    else
        res = 1.0;
    end
end

function move_rotate_rescale(pg, dl, ei, c, T, bta1, bta2, cw1, cw2, E_total)
    if pg.isempty()
        return;
    end
    
    pidx = randsample(1:1:pg.length(),1);
    p = pg.Points(pidx);
    
    while true
        is_move = randsample([true, false],1,true,[1.0/2, 1.0/2]);
        is_rotate = randsample([true, false],1,true,[1.0/2, 1.0/2]);
        if is_move || is_rotate
            break;
        end
    end
    
    pos_new = p.pos;
    ori_new = p.orient;
    
    if is_move
        [pos_new, pro] = pos_noise(p.pos, T+1.0, ei);
        if pro<0
            return;
        end 
    end
    
    if is_rotate
        ori_new = ori_noise(p.orient, (T/p.scale) + 1.0);
    end
    
    p_ = p_copy(p);
    p_.pos = pos_new;
    p_.orient = ori_new;
    
    conn_E_old = p_conn_E(p, bta1, bta2);
    conn_E_new = p_conn_E(p_, bta1, bta2);
    
    E_coll_old = pg.coll_check(p, p, cw1, cw2);
    E_coll_new = pg.coll_check(p_, p, cw1, cw2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bend_E_old = p_bend_E(p, 100.0, 0.5);
    bend_E_new = p_bend_E(p_, 100.0, 0.5);
    
    dE4 = bend_E_new - bend_E_old;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dE3 = E_coll_new - E_coll_old;
    dE1 = vessal_likely(dl, p) - vessal_likely(dl, p_);
    dE2 = conn_E_new-conn_E_old;
    R_move = exp(-(-c*(dE1/T) + dE2/T + dE3/T + dE4/T));
    if ~reject(R_move)
        pg.p_move_rotate_scale(p, pos_new, ori_new);
        
        E_total.E = E_total.E - (E_coll_old+vessal_likely(dl, p_)+conn_E_old+bend_E_old);
        E_total.E = E_total.E + (E_coll_new+vessal_likely(dl, p)+conn_E_new+bend_E_new);
        E_total.E_array(end+1) = E_total.E;
        
        E_total.E_int_array(end+1) = E_total.E_int_array(end) - E_coll_old + E_coll_new - conn_E_old + conn_E_new - bend_E_old + bend_E_new;
        E_total.E_ext_array(end+1) = E_total.E_ext_array(end) - vessal_likely(dl, p_) + vessal_likely(dl, p);
        E_total.p_num_array(end+1) = E_total.p_num_array(end);
    end
end

function track_birth(pg, dl, ei, orient_test, p_scale_min, p_scale_max, bta1, bta2, cw1, cw2, k, c, T, typ, E_total)
    if pg.is_ter_empty()
        return;
    end
    
    pidx = randsample(1:1:size(pg.ter_Points, 2),1);
    p = pg.ter_Points(pidx);
    
    empty_sides = p.head_node_p.get_empty_side();
    
    if size(empty_sides, 2) == 0
        disp("Error in track_birth: Unknow Error!");
        return;
    end
    
    if size(empty_sides, 2) == 2
        side = randsample(empty_sides,1,true,[1.0/2, 1.0/2]);
    else
        side = empty_sides;
    end
    
    x_opt = opt_pos(p, side);
    [x_opt, p_x_opt] = pos_noise(x_opt, (1/2)*(T+0.5*p.scale), ei);
    if p_x_opt<0
        return;
    end
    
    if typ == 1
        ori_new = orient_test(x_opt(2), x_opt(1));
        p_new_ori = 1.0;
    elseif typ == 2
        ori_new = p.orient;
        p_new_ori = 1.0;
    end
    
    if abs(ori_new-p.orient)>=90.0
        side_new = side;
    else
        side_new = -side;
    end
    
    p_birth = point_(x_opt, p_scale_min, p_scale_max);
    p_birth.orient = ori_new;
    p_birth.birth_type = 2;
    
    E_coll = pg.coll_check(p_birth, -1, cw1, cw2);
    E_ext = (vessal_likely(dl, p_birth)-k)*c;
    E_conn = C_Edge(p, side, p_birth, side_new, bta1, bta2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_bend = 100.0*(1.0 - abs(cos(pi*(p.orient-p_birth.orient)/180)) )^(2*0.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dE = (E_coll+E_ext+E_conn+E_bend)/T;
    
    R = exp(-dE) /(p_x_opt*p_new_ori);
    
    if ~reject(R)
        e = node_p(p_birth);
        e.beg = side;
        e.ter = side_new;
        pg.append(p_birth);
        pg.add_conn(p, e, -1);
        
        E_total.E = E_total.E + (E_coll+E_ext+E_conn+E_bend);
        E_total.E_array(end+1) = E_total.E;
        
        E_total.E_int_array(end+1) = E_total.E_int_array(end) + E_coll + E_conn + E_bend;
        E_total.E_ext_array(end+1) = E_total.E_ext_array(end) + E_ext;
        E_total.p_num_array(end+1) = E_total.p_num_array(end) + 1;
    end
end

function track_death(pg, dl, bta1, bta2, cw1, cw2, k, c, T, E_total)
    if pg.is_ter_empty()
        return;
    end
    
    pidx = randsample(1:1:size(pg.ter_Points, 2),1);
    p = pg.ter_Points(pidx);
    
    empty_sides = p.head_node_p.get_empty_side();
    
    if size(empty_sides, 2) == 0
        disp("Error in track_birth: Unknow Error!");
        return;
    end
    
    if size(empty_sides, 2) == 2
        return;
    else
        side = -1*empty_sides;
    end
    
    side_conn = p.head_node_p.get_side_conn(side, true);
    
    if size(side_conn, 2) >1
        disp("Strange in track_death: a side of an end point gets two connections! ");
        return;
    end
    
    p_adj = side_conn.point;
    if p_adj.head_node_p.len<2
        return;
    end
    
    E_coll = pg.coll_check(p, p, cw1, cw2);
    E_ext = (vessal_likely(dl, p)-k)*c;
    E_conn = C_Edge(p, side_conn.beg, p_adj, side_conn.ter, bta1, bta2);
    dE = -(E_coll+E_ext+E_conn)/T;
    
    x_opt = opt_pos(p_adj, side_conn.ter);
    p_x = mvnpdf(p.pos,x_opt,((1/2)*(T+0.5*p.scale))*eye(2));
    
    p_n = 1.0;
    
    R = exp(-dE) * (p_x*p_n);
    
    if ~reject(R)
        pg.dele_conn(p, p_adj);
        pg.dele_(p, pidx);
        
        E_total.E = E_total.E - (E_coll+E_ext+E_conn);
        E_total.E_array(end+1) = E_total.E;
        
        E_total.E_int_array(end+1) = E_total.E_int_array(end) - E_coll - E_conn;
        E_total.E_ext_array(end+1) = E_total.E_ext_array(end) - E_ext;
        E_total.p_num_array(end+1) = E_total.p_num_array(end) - 1;
    end
end

function x_opt = opt_pos(p, side)
    pos = p.pos;
    orient = p.orient;
    scale = p.scale;
    n = [cos((orient/180)*pi), sin((orient/180)*pi)];
    
    x_opt = pos + scale*side*n;
end

function [pos_new, pro] = pos_noise(pos, sig, ei)
    pos_new = floor(mvnrnd(pos,sig*eye(2),1));
    pro = mvnpdf(pos_new,pos,sig*eye(2));
    i = 0;
    while pos_new(2)>size(ei, 1) || pos_new(1)>size(ei, 2) || pos_new(2)<=0 || pos_new(1)<=0
        if i>=10
            pro = -1;
            return;
        end
        pos_new = floor(mvnrnd(pos,sig*eye(2),1));
        pro = mvnpdf(pos_new,pos,sig*eye(2));
        
        i=i+1;
    end
end

function [ori_new, pro] = ori_noise(ori, sig_r)
    n = [cos((ori/180)*pi), sin((ori/180)*pi)];
    n_new = mvnrnd(n,sig_r*eye(2),1);
    
    pro = G_n(n, n_new, sig_r);
    
    ori_new = 180.0*(angle(n_new(1)+1i*n_new(2))/pi);
    if ori_new<0
        ori_new = ori_new+180.0;
    end
end

function res = G_n(n, n_new, sig)
    n1 = n./norm(n);
    n2 = n./norm(n_new);
    cos = dot(n1, n2);
    
    c1 = sig*sqrt(pi/2.0);
    c2 = sig*sqrt(2.0);
    c3 = 2.0*sig^2;
    
    res = c1 * (1.0+erf(cos/c2)) * exp(-(1.0-cos^2)/c3);
end

function res = p_conn_E(p, bta1, bta2)
    res = 0.0;
    nl = p.head_node_p;
    if nl.isempty()
        return;
    end
    
    tem = nl.head;
    while isa(tem,'node_p')
        res = res+C_Edge(p, tem.beg, tem.point, tem.ter, bta1, bta2);
        tem = tem.next_p;
    end
end

function res = p_bend_E(p, c1, c2)
    res = 0.0;
    nl = p.head_node_p;
    if nl.isempty()
        return;
    end
    
    tem = nl.head;
    while isa(tem,'node_p')
        E_bend = c1*(1.0 - abs(cos(pi*(p.orient-tem.point.orient)/180)) )^(2*c2);
        res = res+E_bend;
        tem = tem.next_p;
    end
end

function [col_idx, row_idx] = pos_convert_stack(pos, img_col_len)
    ctest = mod(pos, img_col_len);
    if ctest==0
        row_idx = floor(pos./img_col_len);
        col_idx = img_col_len;
    else
        row_idx = floor(pos./img_col_len)+1;
        col_idx = ctest;
    end
end

function pos_pbirth = pos_convert_flat(pos, img_col_len)
    pos_pbirth = (pos(2)-1)*img_col_len + pos(1);
end

function mid = vessal_likely(dl, p)
    pos = p.pos;
    ori = p.orient;

    low_bound2 = size(dl.theta(dl.theta<=ori), 2);
    
    if low_bound2<=0
        low_bound2 = 1;
    end
    
    if low_bound2>=size(dl.theta, 2) || dl.theta(low_bound2)==ori
        up_bound2 = low_bound2;
    else
        up_bound2 = low_bound2+1;
        ratio2 = (ori-dl.theta(low_bound2))/(dl.theta(up_bound2) - dl.theta(low_bound2));
    end

    aa = dl.AODF_F(low_bound2:up_bound2,pos(2), pos(1));

    if size(aa, 1)>1
        mid = ratio2*aa(1,1) + (1-ratio2)*aa(2,1);
    else
        mid = aa;
    end
end

function res = reject(R)
    if R>=1
        res=false;
        return;
    end
    
    if rand()<= R
        res=false;
    else
        res=true;
    end
end

function res = e_copy(e)
    res = node_p(e.point);
    res.beg = e.beg;
    res.ter = e.ter;
end

function res = p_copy(p)
    res = point_([-1, -1], 0, 1);
    res.pos = p.pos;
    res.orient = p.orient;
    res.scale = p.scale;
    res.head_node_p = p.head_node_p;
end

function res = p_conn(e_id, p_e, e)
    pos = e.point.pos;

    if e_id<0 && (pos(1)>0 || pos(2)>0)
        disp("Index error in connection!");
        res = [];
        return;
    elseif pos(1)>0 || pos(2)>0
        res = p_e(e_id);
        return;
    else
        disp("Unknow error in connection!");
        res = [];
    end
end

function [res, E0, dE] = conn_R(p, e0, e0_idx, e_new, e_new_id, p_e, bta1, bta2, T)
    if e0.point.pos(1)<0 || e0.point.pos(2)<0
        p_e0 = p_e(end);
        E0 = 0.0;
    else
        p_e0 = p_conn(e0_idx, p_e, e0);
        E0 = bta1*(sqrt(dot(p.pos-e0.point.pos, p.pos-e0.point.pos)) - (p.scale+e0.point.scale)/4.0)^2;
    end
    
    p_e_new = p_conn(e_new_id, p_e, e_new);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_new = C_Edge(p, e_new.beg, e_new.point, e_new.ter, bta1, bta2);
    E_bend = 100.0*(1.0 - abs(cos(pi*(p.orient-e_new.point.orient)/180)) )^(2*0.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dE = E_new + E_bend;
    
    res = exp((E0-E_new-E_bend)/T) * (p_e0/p_e_new);
end

function res = C_Edge(p1, side1, p2, side2, bta1, bta2)
    if p2.pos(1)<0 || p2.pos(2)<0
        res = 0.0;
        return;
    end
    
    c1 = C_conn(p1, side1, p2, side2);
    c2 = C_scale(p1, p2);
    res = ((p1.scale+p2.scale)/2.0) * (bta1*c1 + bta2*c2);
end

function res = C_conn(p1, side1, p2, side2)
    x1 = p1.pos;
    x2 = p2.pos;
    orient1 = p1.orient;
    orient2 = p2.orient;
    s1 = p1.scale;
    s2 = p2.scale;
    d1 = s1/4.0;
    d2 = s2/4.0;
    
    x_ = (d2*x1+d1*x2)./(d1+d2);
    
    n1 = d1.*[cos((orient1/180)*pi), sin((orient1/180)*pi)];
    t1 = x1+side1*n1;
    
    n2 = d2.*[cos((orient2/180)*pi), sin((orient2/180)*pi)];
    t2 = x2+side2*n2;
    
    res = dot(x_-t1, x_-t1)/s1^2 + dot(x_-t2, x_-t2)/s2^2;
end

function res = C_scale(p1, p2)
    s1 = p1.scale;
    s2 = p2.scale;
    
    if s1>=s2
        res = (s1/s2)^2-1;
    else
        res = (s2/s1)^2-1;
    end
end

function non_max = non_maxsupress(score, angle_complex, TL, TH)
    [w, h] = size(score);
    
    max_response_ = zeros(size(score));
    
    for r=2:w-1
        for c=2:h-1
            angle_ = angle(angle_complex(r, c))*180/pi;
            
            mag = score(r, c);
            
            if abs(angle_) < 22.5 || abs(angle_) > 157.5
                left = score(r, c - 1);
				right = score(r, c + 1);
                if mag >= left && mag >= right
					max_response_(r, c) = mag;
                end
            end
            
            if (angle_ >= 67.5 && angle_ <= 112.5) || (angle_ >= -112.5 && angle_ <= -67.5)
                top = score(r - 1, c);
				down = score(r + 1, c);
				if mag >= top && mag >= down
					max_response_(r, c) = mag;
                end
            end
            
            if (angle_ > 112.5 && angle_ <= 157.5) || (angle_ > -67.5 && angle_ <= -22.5)
                right_top = score(r - 1, c + 1);
				left_down = score(r + 1, c - 1);
				if mag >= right_top && mag >= left_down
					max_response_(r, c) = mag;
                end
            end
            
            if (angle_ >= 22.5 && angle_ < 67.5) || (angle_ >= -157.5 && angle_ < -112.5)
                left_top = score(r - 1, c - 1);
				right_down = score(r + 1, c + 1);
				if mag >= left_top && mag >= right_down
					max_response_(r, c) = mag;
                end
            end
        end
    end
    
    max_response = dataclass();
    max_response.data = zeros(size(score));
    
    for r=2:w-1
        for c=2:h-1
            mag = max_response_(r, c);
            if mag >= TH
                max_response.data(r, c) = 255;
            elseif mag < TL
                max_response.data(r, c) = 0;
            end
        end
    end
    
    non_max = max_response.data;
    
end