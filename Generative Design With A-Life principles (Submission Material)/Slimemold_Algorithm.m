close all
F = findall(0,'type','figure','tag','TMWWaitbar')
delete(F)

[welcome_message,ploters,video,popup,video_name,multi_put] = start_prompt()
welcome(welcome_message)

[model,coords,lattice,results,COM_index,COM_loc,dis_ma,t] = import_gen_model('.\3D Geometry\3DBracket.stl',ploters)

seed = COM_index % put the seed for the slime mold in the centre of mass for the part, this is a continued practice from the swarm method

angle1_Cell  = {0} % set up all angles and positions so that the seed can then be iterated upon
angle2_Cell  = {0}

record = []
cell = {seed}
counter = 1

glob_loc_pos  = seed % each iteration has a local list of positions and angles which is added to the global lists which hold each itereative list of variabls in different cells
glob_loc_ang  = 0
glob_loc_ang2 = 0

if video == 1 % set up video plotting if selected
    figure  

    if exist('video_num','var') == 1
        video_num = video_num + 1

    else
        video_num = 0
    end

    name_vid = string(video_num)

    myVideo = VideoWriter(name_vid,'MPEG-4'); 
    myVideo.FrameRate = 10; 
    open(myVideo)
end

iters = 200 

% Nested loops to iterate on each branch, every iteration, the previous endpoints are given their branch designations for 
% successive branch number and angle
for i = 1:iters 
    
    angle_dump = []
    local_pois  = []
    local_angs  = []
    local_angs2 = []
    
    prev_cell  = cell{:,end}
    prev_angl  = angle1_Cell{:,end};
    prev_angl2 = angle2_Cell{:,end};
    

    for j = 1:length(cell{:,end})

        num_bran        = num_branch(prev_cell(j),lattice,i,iters);
        

        for k = 1:num_bran

            
            dists           = dist_func(prev_cell(j),lattice,model);
            [angle1,angle2] = angle_func(prev_cell(j),lattice,i);
            [new_poi,new_ang1,new_ang2]        = gen_points(prev_cell(j),dists,angle1,angle2,coords,prev_angl(j),prev_angl2(j));
            
            local_pois  = [local_pois,new_poi];
            local_angs  = [local_angs,new_ang1];
            local_angs2 = [local_angs2,new_ang2];
            
            angle_dump = [angle_dump,angle1]

            
            if video == 1
            
                prev_points = coords(prev_cell(j),:)
                prev_points(:,1)

                points_fin = coords(new_poi,:)

                plot3([prev_points(:,1) points_fin(:,1)],[prev_points(:,2) points_fin(:,2)],[prev_points(:,3) points_fin(:,3)])
                hold on
                view(45,45)

                drawnow 

                pause(0.01) 

                frame = getframe(gcf); %get frame
                writeVideo(myVideo, frame);
                 
            end

        end

    end
    
    cell{1,i+1}        = local_pois
    angle1_Cell{1,i+1} = local_angs
    angle2_Cell{1,i+1} = local_angs2


    counter = counter + 1  

end

if video == 1
    close(myVideo)
end

index_Arr = cell2mat(cell)
points_final = coords(index_Arr,:) % retrieve slime mold activated nodes from the part coordinate array
figure
scatter3(points_final(:,1),points_final(:,2),points_final(:,3),10,'r') % scatter plot to show final structure

save('slime_mold_results.mat','points_final') % Export the results as a matrix file for the STL exporter file to convert into an STL outfile

function [welcome_message,ploters,video,popup,video_name,multi_put] = start_prompt()

prompt = {'Would you like to plot video) ?','Would you like to see plots ?','Would you like a welcome message ?','If you would like a custom video made, please enter name here' ,'Did you want a multi_ouputsystem ?'};
dlgtitle = 'Input, please enter 1 in boxes you wish to see';
dims = [1 35];
definput = {'0','0','0','0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);


welcome_message = str2num(answer{3});
ploters = str2num(answer{2});
video   = str2num(answer{1});
popup   = 1;
video_name = answer{4};
multi_put = str2num(answer{5})

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
function welcome(y_n)
if y_n == 1
    msgbox('The slime mold is based on physarum polycephalum, a monocellular organism that can expand out and find a food source, here using forward kinematics, its tendrils or branches can be simulateded expanding out into the nodes of the part. If you have changed the iterations from 75 to a higher number please be aware that the time for completion ( if video plotting is selected ) will be high.')
end
end % Display welcome message (if user wants) 

function centroid = find_cent(points)

    centroid = [mean(points(:,1)),mean(points(:,2)),mean(points(:,3))];
  
end % find centre of mass of part
function [k,dist] = cit_group(citizen,group)
[k,dist] = dsearchn(group,citizen);
end % Find closest part node to slime tendril which will become new position for new endpoint



% Here we take the lattice input (nutritional value) described in the
% technical report, this lattice input with then, based on measured slime
% mold values, decide the branch number, distance of tendril segment, angle
% of bifurication 

function branches        = num_branch(input,lat,i,iters)

if i ==1 
    branches = 50
    
else
input = lat(input)

p0 = (0.5/(1 + exp((0.4*20*input)-(0.4*10)))         )*100
p1 = (0.2/(1 + exp((0.5*20*input)-(0.5*10))) + 0.25  )*100
p2 = (0.6/(1 + exp(-(0.4*20*input)-(0.4*10))) + 0.04 )*100
p3 = 100 - (p0+p1+p2)                                

p0 = round(p0)
p1 = round(p1)
p2 = round(p2)
p3 = round(p3)

prob_arr = horzcat(ones(1,p0)*1, ones(1,p1)*1 , ones(1,p2)*1  , ones(1,p3)*2)

marker = randi([1 length(prob_arr)])

branches = prob_arr(marker)
end

end                                                   
function dis             = dist_func(input,lat,model)

min_d = model.Mesh.MinElementSize
max_d = model.Mesh.MaxElementSize
input = lat(input)

dis = 0.25*(max_d-min_d)*(log((0.5*input)/(10-10*input))) + max_d

if dis >=1% change to find max-min gemoetry values
    dis = 1% change to find max-min gemoetry values
elseif dis <= 0.1 % change to find max-min gemoetry values
    dis = 0.1% change to find max-min gemoetry values
end
end % Here the distance of a tendril (branch) segment is, based on measures slime mold variables  
function [an1,an2]       = angle_func(input,lat,iter)

if iter == 1
    
    an1 = randi([-360,360])
    an2 = randi([-360,360])
    
else
    input = lat(input)

    lower = 10/(1+exp(-((10*input)-(0.5*10)))) % 45 changed to 10
    upper = 20/(1+exp(-((10*input)-(0.5*10)))) + 10 % 80>20, 45 > 10

    an1 = randi([round(lower) round(upper)])
    an2 = randi([round(lower) round(upper)])
end

an1 = (2*pi*an1)*neg_pos()/360
an2 = (2*pi*an2)*neg_pos()/360

end    
function [n_p,n_a1,n_a2] = gen_points(p_c,di,an1,an2,cords,prev_angl,prev_angl2)

old_cord = cords(p_c,:)

n_a1 = an1 + prev_angl
n_a2 = an2 + prev_angl2

x_n = old_cord(1) + cos(n_a1)*di
y_n = old_cord(2) + sin(n_a1)*di
z_n = old_cord(3) + sin(n_a2)*di

n_p_appx = [x_n,y_n,z_n]

n_p = knnsearch(cords,n_p_appx)



end
function nplus= neg_pos()
random = randi([1,100])

if random > 50
    nplus = 1
else
    nplus = -1
end
end