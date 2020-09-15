
prompt = 'Please enter either > slime_mold < or > swarm_optimised < depending on the result types '
in_out_name = input(prompt','s')

loadm = strcat(in_out_name,'_results.mat')

load(loadm)

[X,Y,Z] = sphere(5);

X = X(:) ;
Y = Y(:) ;
Z = Z(:) ;

scatter3(X,Y,Z)

[model,p,e,t] = import_gen_model('.\3D Geometry\3DBracket.stl')
scatter3(points_final(:,1),points_final(:,2),points_final(:,3))


sph_tot = []

for i = 1:length(points_final)
    
    x_T = X + points_final(i,1);
    y_T = Y + points_final(i,2);
    z_T = Z + points_final(i,3);
    
    disp(i)
    sph_tot = [sph_tot;[x_T,y_T,z_T]];
    
    
end

[k,dist] = dsearchn(sph_tot,p)

idx_move = find(dist > 0)
k_idx_move = k(idx_move)
p(idx_move,:) = sph_tot(k_idx_move,:)

figure
scatter3(p(:,1),p(:,2),p(:,3),'r')

facets = [t(:,[1 2 3]);t(:,[1 2 4]);t(:,[1 3 4]);t(:,[2 3 4])];
facets = sort(facets,2);
facets = sortrows(facets);
duploc = find(all(diff(facets) == 0,2));
facets([duploc;duploc + 1],:) = [];

TR = triangulation(facets,p)

savep = strcat(in_out_name,'_part.stl')

stlwrite(TR,savep)

function [model,p,e,t] = import_gen_model(name)

    model = createpde('structural','static-solid');
    importGeometry(model,name);
    generateMesh(model);
        
    figure
    pdegplot(model,'FaceLabels','on')
    view(-134,-32)
    title('Bracket with Face Labels, Rear View')
    
    [p,e,t] = meshToPet(model.Mesh);
    p = transpose(p);
    t = transpose(t);

end 