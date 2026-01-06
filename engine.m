clc
close all
clearvars



% Input dimension
Lx = inputDefault('Domain Size [Lx = Ly]',100);
Ly = Lx;
Nx = inputDefault('Grid Resolution [Nx = Ny]',101);
Ny = Nx;
epsilon = inputDefault('Permittivity',1);

% Initialize Potential & Charge Spaces
h = Lx/(Nx - 1);
V = zeros(Nx, Ny);
rho = zeros(Nx, Ny);

% Interactive Charge Editor
figure('Name','Click to Place Charges','NumberTitle','off');
imagesc(V'); axis equal tight; colorbar;axis xy;
title('Click to place charges (right click to finish)');
xlabel('X'); ylabel('Y');
hold on;

charge_done = false;
while ~charge_done
    [x, y, button] = ginput(1);  % get one click
    if button == 3   % right-click to finish
        charge_done = true;
        break;
    end
    q = input('Enter charge value at this point: ');

    % Convert clicked coordinates to grid indices
    i = round(x); % X maps to first index (Rows in V, X in V')
    j = round(y); % Y maps to second index (Cols in V, Y in V')
    
    if i >= 1 && i <= Nx && j >=1 && j <= Ny
        rho(i,j) = rho(i,j) + q / (h^2); % Equation 10
        plot(i, j, 'ro', 'MarkerSize',10, 'LineWidth',2); 
        text(i, j, num2str(q), 'Color','w','FontWeight','bold'); 
    else
        disp('Clicked outside grid.');
    end
end


% Boundary Condition
% RIGHT wall
RightWallMode = inputDefault('Right Wall Dirichlet / Neumann [1/2]',1);
if RightWallMode == 1
    RightWallPotential = inputDefault('Potential',0);
    V(:,Ny) = RightWallPotential;
else
    RightWallGn = inputDefault('Gradient',0);
end

% LEFT wall
LeftWallMode = inputDefault('Left Wall Dirichlet / Neumann [1/2]',1);
if LeftWallMode == 1
    LeftWallPotential = inputDefault('Potential',0);
    V(:,1) = LeftWallPotential;
else
    LeftWallGn = inputDefault('Gradient',0);
end

% TOP wall
TopWallMode = inputDefault('Top Wall Dirichlet / Neumann [1/2]',1);
if TopWallMode == 1
    TopWallPotential = inputDefault('Potential',0);
    V(1,:) = TopWallPotential;
else
    TopWallGn = inputDefault('Gradient',0);
end

% BOTTOM wall
BottomWallMode = inputDefault('Bottom Wall Dirichlet / Neumann [1/2]',1);
if BottomWallMode == 1
    BottomWallPotential = inputDefault('Potential',0);
    V(Nx,:) = BottomWallPotential;
else
    BottomWallGn = inputDefault('Gradient',0);
end




% Method, Tolerance & Maximum Iteration
method = inputDefault('Jacobi / Gauss–Seidel / SOR [1/2/3]', 3);
if method == 3
    omega = inputDefault('w for SOR [1.5 ~ 1.9]', 1.8);
end
tolerance =inputDefault('Tolerance', 1e-6);
maxIter = inputDefault('Maximum iteration',20000);
errorTrack = zeros(maxIter,1);


% GUI Figure
figure('Name','2D Poisson Solver','NumberTitle','off');
subplot(1,2,1);
hPlot = imagesc(V'); axis equal tight; colorbar;axis xy;
title('Potential Distribution \phi')
xlabel('X'); ylabel('Y');

subplot(1,2,2);
hold on;
hErr = plot(errorTrack);
xlabel('Iteration'); ylabel('Error');
title('Convergence Error');
ylim([10^-7  10]);
set(gca,'YScale','log');


% Solve

for iter = 1:maxIter
    % Neuman Update
    % LEFT wall
    if LeftWallMode == 2
        V(:,1) = V(:,2) + h * LeftWallGn;
    end

    % RIGHT wall
    if RightWallMode == 2
        V(:,Ny) = V(:,Ny-1) + h * RightWallGn;
    end

    % TOP wall
    if TopWallMode == 2
        V(1,:) = V(2,:) + h * TopWallGn;
    end

    % BOTTOM wall
    if BottomWallMode == 2
        V(Nx,:) = V(Nx-1,:) + h * BottomWallGn;
    end


    % Jacobian
    if method == 1
        V_temp = V;
        V_old = V;
        for i = 2:Nx-1
            for j = 2:Ny-1
                V_temp(i,j) = 0.25 * (V_old(i+1,j) + V_old(i-1,j) + V_old(i,j+1) + V_old(i,j-1) + (h^2/epsilon)*rho(i,j));
            end
        end
        V = V_temp;

    % Gauss-Siedal
    elseif method == 2
        V_old = V;
        for i = 2:Nx-1
            for j = 2:Ny-1
                V(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1) + (h^2/epsilon)*rho(i,j));
            end
        end

    % SOR
    else
        V_old = V;
        for i = 2:Nx-1
            for j = 2:Ny-1
                V_star = 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1) + (h^2/epsilon)*rho(i,j));
                V(i,j) = (1-omega)*V(i,j) + omega*V_star;
            end
        end
    end



    % Re-impose Dirichlet boundary nodes

    % RIGHT wall
    if RightWallMode == 1
        V(:,Ny) = RightWallPotential;
    end

    % LEFT wall
    if LeftWallMode == 1
        V(:,1) = LeftWallPotential;
    end

    % TOP wall
    if TopWallMode == 1
        V(1,:) = TopWallPotential;
    end

    % BOTTOM wall
    if BottomWallMode == 1
        V(Nx,:) = BottomWallPotential;
    end


    % Track error
    diff = max(max(abs(V - V_old)));
    errorTrack(iter) = diff;

    % Update GUI every 20 iterations
    if mod(iter,20) == 0
        set(hPlot,'CData',V');
        set(hErr,'YData',errorTrack);
        drawnow;
    end

    % Convergence check
    if diff < tolerance
        fprintf("Converged in %d iterations\n",iter);
        break
    end
end

% Electric Field Map
[Ex, Ey] = gradient(-V, h);
figure('Name','Electric Field','NumberTitle','off');
quiver(Ey',Ex');axis xy
axis equal tight;
title('Electric Field Vectors');
xlabel('X'); ylabel('Y');


% Potential Map
figure('Name','Final Potential','NumberTitle','off');
surf(V'); shading interp; colorbar;axis xy
title('3D Potential Distribution');
xlabel('X'); ylabel('Y'); zlabel('\phi');


% Helper Function for default input
function val = inputDefault(prompt, default)
str = input([prompt ' (default ' num2str(default) '): '], 's');
if isempty(str)
    val = default;
else
    val = str2double(str);
end
end


