close all
F = findall(0,'type','figure','tag','TMWWaitbar')
delete(F)

[welcome_message,ploters,video,popup,video_name] = start_prompt()

welcome(welcome_message)

[model,p,lattice,results,COM_index,COM_loc,dis_ma,t] = import_gen_model('.\3D Geometry\3DBracket.stl',ploters)


new_pop = seed_earth(COM_index,p,0.2);
size_pop = size(new_pop)
size_pop = size_pop(1)
movement_matrix = zeros(size_pop,1)

if ploters == 1
    plotter(p,COM_loc,'Complete set of mesh points (with centre)','b')
end

ite = 10000;

tic
result_swarm = iteration_function(ite,ploters,video,popup,new_pop,p,lattice,movement_matrix,video_name)
toc


points_final = result_swarm
save('Swarm_optimised_results.mat','points_final')

% Introductory functions, setup, import and ready up geometry

function [welcome_message,ploters,video,popup,video_name] = start_prompt()

prompt = {'Would you like to plot video) ?','Would you like to see plots ?','Would you like a welcome message ?','If you would like a custom video made, please enter name here' };
dlgtitle = 'Input, please enter 1 in boxes you wish to see';
dims = [1 35];
definput = {'0','0','0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);


welcome_message = str2num(answer{3});
ploters = str2num(answer{2});
video   = str2num(answer{1});
popup   = 1;
video_name = answer{4};

end % Get user input on running preferences and UI
function [model,p,lattice,result,COM_index,COM,Dist_matrix,tetra] = import_gen_model(name,ploters)

    model = createpde('structural','static-solid');
    importGeometry(model,name);
    generateMesh(model);
    structuralProperties(model,'YoungsModulus',200e9, ...
                               'PoissonsRatio',0.3)            
    structuralBC(model,'Face',12,'Constraint','fixed');

    structuralBoundaryLoad (model,'Face',17,'SurfaceTraction',[0;0;-1E11]);
    result = solve(model);
    
    if ploters == 1
    
        figure
        pdegplot(model,'FaceLabels','on')
        view(-134,-32)
        title('Bracket with Face Labels, Rear View')

        figure
        pdeplot3D(model,'ColorMapData',result.VonMisesStress)
        title('x-displacement')
        colormap('jet')
    end

    vonMis = result.VonMisesStress;
    scaled_inverted_vM = 1 - rescale(vonMis);

    [p,e,t] = meshToPet(model.Mesh);
    p = transpose(p);
    tetra = transpose(t)

    lattice = zeros(size(scaled_inverted_vM));
    lattice = lattice + scaled_inverted_vM ;
    
    COM = find_cent(p);
    [COM_index,COM_dist]    = cit_group(COM ,p); 
    [All_index,Dist_matrix] = cit_group(p, COM);
    
    Dist_matrix = 1 - rescale(Dist_matrix);
    
    lattice = lattice .* Dist_matrix;

end % Import the model, gen mesh and give PDE results
function f = loading_bar(it_loc,it_glob,f)
    if it_loc == 1;
        f = waitbar(it_loc/it_glob,'Please wait');
    end
    
    if ismembertol(it_loc/it_glob, [linspace(0.01,1,100)],1E-5);
        close(f);
        f = waitbar(it_loc/it_glob,sprintf('Please wait, %d %% Complete',100*it_loc/it_glob));
    end
end % Provide a loading bar to give progress feedback
function welcome(y_n)
if y_n == 1
    msgbox('Good day, this is a summary of the variables in the code so far .The general structure is thusly: We have four structure types ,The lattice, the swarm, the inhabitant matrix and the mesh nodes. The lattice is the available areas the swarm may travel to, this means that is influenced by stress values at that point . The swarm is also influenced by the inhabitance matrix, which keeps an up to date record of where swarm members are, ensuring that no node is double parked by two swarm members. The swarm are the number of nodes that will represent our new structure, over time they will travel to areas of high stress and low swarm density to not only fill the shape but emphasise material placement where most needed. The stresses are Von misses, stress values corresponding to a mesh that represents the geometry of the part.')
end
end % Display welcome message (if user wants) 


% Main looping function, runs through and moves swarm members over a
% specified number of iterations
function return_swarm_points = iteration_function(iters,ploters,video,popup,new_pop,p,lattice,movemat,vid_name)

figure

if video == 1
    
    if vid_name == 0
        nom = string(datetime('today'))
    else
        nom = vid_name
    end
    
    myVideo = VideoWriter(nom,'MPEG-4'); %open video file
    myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    open(myVideo)
    
end

f = 0;

for i = 1:iters

    f = loading_bar(i,iters,f);
    
    [index_place,dist] = cit_group(new_pop,p);
    
    pop_prese = index_place;
    
    [minval,new_index] = min(movemat);
    movemat(new_index) = movemat(new_index) + 1;

    
    [new_pop,lattice] = lookabout(lattice,new_index,new_pop,p,pop_prese);
    
    if video == 1
        
        if rem(i,100) == 0
            scatter3(new_pop(:,1),new_pop(:,2),new_pop(:,3),10,'b')
            view(i/100,45)

            xlim([min(p(:,1)) max(p(:,1))])
            ylim([min(p(:,2)) max(p(:,2))])
            zlim([min(p(:,3)) max(p(:,3))])

            drawnow 

            pause(0.01) 

            frame = getframe(gcf); %get frame
            writeVideo(myVideo, frame);
        end
    end
end
if video == 1
    close(myVideo)
end

if ploters == 1
    cen = find_cent(new_pop)
    plotter(new_pop,cen,sprintf('Colony or swarm after %d iterations (with centre)',iters),'y')
end
return_swarm_points = new_pop

end

% Querey/Action functions, checks lattices and makes movements based on
% criteria
function [population_Return,lat_return] = lookabout(lattice,ind,new_pops,coords,pop_pres)

population_Return = new_pops;

MinValPrev = 0;

present = 1;
[individual,distance] = cit_group(new_pops(ind,:),coords);

counter = 0;

while present == 1 
    
    num_neighbours = 9;
    
    if counter > 2 & MinIndex == MinValPrev;
        lattice(MinIndex) = lattice(MinIndex) + 1;
        
    elseif counter > 500
        
        num_neighbours = 20;
    end  
    
    
    
    [mIdx,mD] = knnsearch(coords,coords(individual,:),'k',num_neighbours);
    
    mIdx = mIdx(2:9);
    
    vals = lattice(mIdx,:);

    [MinValue,MinIn] = min(vals);
    
    MinIndex = mIdx(MinIn);
    
    [individual,distance] = cit_group(coords(MinIndex,:),coords);
    
    counter = counter + 1;
    
    if ~ismember(MinIndex,pop_pres);
        
        present = 0;
    end

    

    
    MinValPrev = MinIndex;
   
end

lat_return = lattice;
population_Return = move(ind,new_pops,coords,MinIndex);

end % Primary source of behaviour, looks for empty lattice slots and picks lowest values
function plotter(object,center,titl,col)

figure
scatter3(object(:,1),object(:,2),object(:,3),10,col)

if center ~= 0
    hold on
    scatter3(center(:,1),center(:,2),center(:,3),1000,'*g')

end
title(titl)
end % Plots swarm as series of points and provides graph 
function [k,dist] = cit_group(citizen,group)
[k,dist] = dsearchn(group,citizen);
end % Here the citezen of the swarm checks its closest connection to group, to keep a cohesive structure
function population = seed_earth(seed,coords,required_density)
    
    num_neighbours = round(required_density*length(coords));

    [mIdx,mD] = knnsearch(coords,coords(seed,:),'k',num_neighbours);
    
    population = coords(mIdx,:);
    
end % Here we poll the coordinates to see the landscape with swarm population
function centroid = find_cent(points)

    centroid = [mean(points(:,1)),mean(points(:,2)),mean(points(:,3))];
  
end % to track how the swarm is moving we will monitor the centre of mass, and try and compare it to the coordinate COM
function return_swarm = move(swarm_index,swarm,coords,new_pos_index)% activate/occupy the nodes of the new swarm location matrix
    
    return_swarm = swarm;
    return_swarm(swarm_index,:) = coords(new_pos_index,:);
    
end

