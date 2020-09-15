h = figure;

xlim([0 6])
ylim([0 6])

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

auxim_list= rand(1500,2)*5
vein_nodes=[[0.1 0.1]];

scatter(vein_nodes(:,1),vein_nodes(:,2),5,'b','filled')
hold on
scatter(auxim_list(:,1),auxim_list(:,2),2,'r','filled')
hold on

for j = 1:60
    
[vein_index,dist_vein]  = dsearchn(vein_nodes,auxim_list)

[b,m1,n1] = unique(vein_index,'first');
[c1,d1] =sort(m1);
b = b(d1);

disp('vein_index,vein_repeat')
disp(vein_index)
disp(b)

new_nodes = [];


    for i = 1:length(b)

        D = 0.1;

        vein_pointer = b(i)

        X = find(vein_index == vein_pointer);
        local_auxim = auxim_list(X,:);

        n_vect = sum((local_auxim - vein_nodes(vein_pointer,:))/(norm(local_auxim - vein_nodes(vein_pointer,:))),1)

        disp('vein pointypoint')
        disp(vein_nodes(vein_pointer,:))

        new_node = vein_nodes(vein_pointer,:) + D*(n_vect/norm(n_vect));

        new_nodes = [new_nodes;new_node];
        
        if j == 1
            disp('pass')
        else
        
            plotter_point = vein_nodes(vein_pointer,:)
            plot([plotter_point(:,1),new_node(:,1)],[plotter_point(:,2),new_node(:,2)])
       

        end
        
    end
    
                drawnow 
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            outfile = 'leaf_pattern.gif';

            % On the first loop, create the file. In subsequent loops, append.
            if j==1
                imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
            else
                imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
            end
        
    vein_nodes = [vein_nodes;new_nodes];
        
    end

      

