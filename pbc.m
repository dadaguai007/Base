%PBC
function E = pbc(Ex,Ey,theta)
%Not tested ， the function is done
   if nargin < 2
        theta = 0;  % Default rotation angle
    end
    
    % Check and reshape E if necessary
    if size(Ex, 2) == 1
        Ex = [Ex, zeros(size(Ex))];  % Assume Y component is 0
    end
        % Check and reshape E if necessary
    if size(Ey, 2) == 1
        Ey = [zeros(size(Ey)),Ey];  % Assume Y component is 0
    end
    
    % Create the rotation matrix
    rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    % inv
    rot_inv1 = inv(rot);
%     rot_inv = rot.';
    rot_inv = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    
    % N × 2 to 2×N
    E_x = [Ex(:,1),Ey(:,1)].';
    E_y = [Ex(:,2),Ey(:,2)].';
    % Rotate the input optical field
    % 1*N 
    Ex_out = rot_inv(1,:)*E_x;
    Ey_out = rot_inv(2,:)*E_y;
    


    % Extract X and Y polarization components
    % 2×N to N × 2
    E = [Ex_out,Ey_out].';
end