
[model,p,e,t] = import_gen_model('.\3D Geometry\FlatPlate1.stl',(0.5));

thinning_index = find(p(:,2) > 0)
p_thin_select = p(thinning_index,:)
p_thin_select(:,2) = p_thin_select(:,2) * 0.75
p(thinning_index,:) = p_thin_select

gx_max = max(p(:,1))
gy_max = max(p(:,3))

y_min   = min(p(:,3))

scale_const = max([gx_max,gy_max-y_min])/5


x_max  = 0.5*scale_const
y_max  = min(p(:,3)) + 0.5*scale_const


stop_flag = 0

vein_nodes = [[0,min(p(:,3))];[0,min(p(:,3))+0.01];[0,min(p(:,3))+0.02]]
auxin_list = [[0,0]]

threshold = 0.02*scale_const
auxin_list_test = auxin_gen(auxin_list,vein_nodes,threshold,y_min,x_max,y_max,5);

counter = 0

figure

while stop_flag == 0 
    
    auxin_list = auxin_gen(auxin_list,vein_nodes,threshold,y_min,x_max,y_max,5)
    
    vein_nodes = nvect_gen(auxin_list,vein_nodes,scale_const);
    
    auxin_list = auxin_kill(auxin_list,vein_nodes,threshold)
    
    [x_max,y_max,stop_flag] = area_expansion(gx_max,gy_max,vein_nodes,x_max,y_max,y_min,scale_const);
    
    counter = counter + 1
    
    if counter == 1000
        stop_flag = 1
    end

end

scatter(vein_nodes(:,1),vein_nodes(:,2),2,'b','filled')
scatter(auxin_list(:,1),auxin_list(:,2),2,'r','filled')

part_alteration(vein_nodes,p,t)

function part_alteration(vein_nodes,p,t)
surf_vert = 0

vein3d = [vein_nodes(:,1),(ones(length(vein_nodes),1)*surf_vert),vein_nodes(:,2)]


idx_surf = find(1e-04 > p(:,2) >= 0)
p_surf = p(idx_surf,:)

figure
scatter3(p_surf(:,1),p_surf(:,2),p_surf(:,3))

k = knnsearch(p_surf,vein3d,'k',5)
k = k(:)

selected = p_surf(k,:)
selected(:,2) = selected(:,2) - 1
p_surf(k,:) = selected

p(idx_surf,:) = p_surf

figure
scatter3(p(:,1),p(:,2),p(:,3))

facets = [t(:,[1 2 3]);t(:,[1 2 4]);t(:,[1 3 4]);t(:,[2 3 4])];
facets = sort(facets,2);
facets = sortrows(facets);
duploc = find(all(diff(facets) == 0,2));
facets([duploc;duploc + 1],:) = [];

TR = triangulation(facets,p)
stlwrite(TR,'leaf_optimised_part.stl')
end % take nodes and warp resulting geometry then export it as STL file

function [model,p,e,t] = import_gen_model(name,mesh_siz)

    model = createpde('structural','static-solid');
    importGeometry(model,name);
    generateMesh(model,'Hmax',mesh_siz);
        
    figure
    pdegplot(model,'FaceLabels','on')
    view(-134,-32)
    title('Bracket with Face Labels, Rear View')
    
    [p,e,t] = meshToPet(model.Mesh);
    p = transpose(p);
    t = transpose(t)

end % Import the model, gen mesh and give PDE results


% As the leaf veins approach the boarders of the current area, if they exceed 75% of teh height or length of these borders, 
% this function expands them out, however, if the borders exceed the maximum dimensions of the part they will no longer be expanded
function [x_ma,y_ma,stop_flag] = area_expansion(glob_x_ma,glob_y_ma,vein_nodes,curre_x,curre_y,curre_ymin,scale_const)

if any(vein_nodes(:,1) > (curre_x*0.75)) | any(vein_nodes(:,1) < (-curre_x*0.75))
    
    if curre_x >= glob_x_ma
        curre_x = glob_x_ma
    else
        curre_x = curre_x + 0.2*scale_const
    end
    
end
    
if any(abs(vein_nodes(:,2))>abs(curre_y*0.75))
    
    if curre_y >= glob_y_ma
        curre_y = glob_y_ma
    else
        curre_y = curre_y + 0.2*scale_const

    end
    
end

if max(vein_nodes(:,2)) >= glob_y_ma & max(vein_nodes(:,1)) >= glob_x_ma
   
    stop_flag = 1
else
    stop_flag = 0
end

plotter(curre_x,curre_ymin,curre_y)

x_ma = curre_x
y_ma = curre_y

end

% Here the closest auxin particles are aissgned to their nearest vein node
% neighbour and the sum of the euclidean norms are calcualted to find the
% new position of a vein node
function vein_nodes_return = nvect_gen(auxin_list,vein_nodes,scale_const)
[vein_index,dist_vein]  = dsearchn(vein_nodes,auxin_list)

[b,m1,n1] = unique(vein_index,'first');
[c1,d1] =sort(m1);
b = b(d1);

disp('vein_index,vein_repeat')
disp(vein_index)
disp(b)

new_nodes = [];

for i = 1:length(b)
    
    D = 0.01*scale_const;
    
    vein_pointer = b(i)
    
    X = find(vein_index == vein_pointer);
    local_auxin = auxin_list(X,:);
    
    n_vect = sum((local_auxin - vein_nodes(vein_pointer,:))/(norm(local_auxin - vein_nodes(vein_pointer,:))),1)
        
    disp('vein pointypoint')
    disp(vein_nodes(vein_pointer,:))
    
    new_node = vein_nodes(vein_pointer,:) + D*(n_vect/norm(n_vect));
    
    new_nodes = [new_nodes;new_node];
end

vein_nodes = [vein_nodes;new_nodes];

vein_nodes_return = vein_nodes

end

% here, as the vein nodes have expanded, the auxin nodes in teh vicinity of
% the leaf vein nodes need to be removed 
function auxin_list_return = auxin_kill(auxin_list,vein_nodes,threshold)

    [k,dist_vein]  = dsearchn(vein_nodes,auxin_list);

    distindex = find(dist_vein < threshold);

    auxin_list(distindex,:) = []
    
    auxin_list_return = auxin_list
end

% Here based on a minimum distance from existing auxin and leaf vein nodes
% new auxin particles are generted within the bounds as defined by the area
% expansion function
function auxin_list_return = auxin_gen(auxin_list,vein_nodes,threshold,y_min,x_max,y_max,new_nodes)
    
    for i = 1:new_nodes
        
        n_x = unifrnd(-x_max,x_max)
        n_y = unifrnd(y_min,y_max)
        
        new_point = [n_x,n_y]
        disp('new_point')
        disp(new_point)
        
        [k1,dist_auxin] = dsearchn(auxin_list,new_point);
        [k2,dist_vein]  = dsearchn(vein_nodes,new_point);
        
        if dist_auxin > threshold & dist_vein > threshold 
            
            auxin_list = [auxin_list; new_point]

        end
 
    end
    
    auxin_list_return = auxin_list
    
end

% plot the resulting data
function plotter(x_max,y_min,y_max)

plot([-x_max,-x_max],[y_min,y_max])
hold on
plot([x_max,x_max],[y_min,y_max])
hold on

plot([-x_max,x_max],[y_min,y_min])
hold on
plot([-x_max,x_max],[y_max,y_max])

end