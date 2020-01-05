function [Usr_distances, Usr_angles, l_array] = Cell_User_deployment(R, N)
%%
%% R: cell size
%% N: number of users to be deployed
%% Output: distances of each user deployed (randomly) within the cell from the BS
%%
figure('Name','Cell Deploy');
rho = zeros(1,N);
theta = zeros(1,N);
for i=1:N
    rho(i) = R*rand();
    theta(i) = 2*pi*rand();
    
    [x,y] = pol2cart(theta(i),rho(i));
    scatter(x,y,'ko');
    text(x,y,sprintf(' %d',i))
    hold on
end
circle(0,0,R);
hold on
plot(0,0,'g+');
xlabel('R');
ylabel('R');
title('User Allocation in the Cell');
grid

Usr_distances = rho;
Usr_angles = theta;

l_array = [];
for i=1:N
    for j = 1:N
        if (i ~= j) && (ismember(i,l_array) == 0) && (ismember(j,l_array) == 0)
            if ((Usr_angles(i) + pi/12) >= Usr_angles(j)) &&  ((Usr_angles(i) - pi/12) <= Usr_angles(j))
                disp('Group paired');
                X = sprintf('User %d with distance: %d. and angle: %d',i,Usr_distances(i), Usr_angles(i));
                fprintf(X);              
                [x,y] = pol2cart(Usr_angles(i),Usr_distances(i));
                usr_i = sprintf('\nUser %d coordinates: ',i);
                fprintf(usr_i);
                disp([x,y]);
                scatter(x,y,'ro');
                               
                Y = sprintf('User %d with distance: %d. and angle: %d',j,Usr_distances(j), Usr_angles(j));
                fprintf(Y);
                [x,y] = pol2cart(Usr_angles(j),Usr_distances(j));
                scatter(x,y,'ro');
                usr_j = sprintf('\nUser %d coordinates: ',j);
                fprintf(usr_j);
                disp([x,y]);
                
                l_array = [l_array i j];
                disp(l_array);
            end
        end
    end
end
hold off

