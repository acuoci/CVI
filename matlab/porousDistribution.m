function [bin_x, bin_y] = porousDistribution(planar, nx,x, ny,y, phi)

max_radius = 9;
n_bins = 36;
delta_bin = max_radius/(n_bins-1);

bins = 0:delta_bin:max_radius;
sum_bins = zeros(n_bins-1,1);

for i=2:1:nx-1
    
    ri = x(i) - 0.5*(x(i)-x(i-1));
    re = x(i) + 0.5*(x(i+1)-x(i));
    
    if (planar == false)
        A = pi*(re^2-ri^2);
    else
        A = re-ri;
    end
    
    for j=2:1:ny-1
        
        hm = y(j) - 0.5*(y(j)-y(j-1));
        hp = y(j) + 0.5*(y(j+1)-y(j));
        
        V = (hp-hm)*A;
        
        for k=2:1:n_bins+1
            if (phi(i,j) < bins(k))
                sum_bins(k-1,1) = sum_bins(k-1,1)+V;
                break;
            end
        end
        
    end
    
    %North
    deltay = 0.50*(y(ny)-y(ny-1));
    V = A*deltay;
    
    for k=2:1:n_bins+1
        if (phi(i,ny) < bins(k))
            sum_bins(k-1,1) = sum_bins(k-1,1)+V;
            break;
        end
    end
    
    %South
    deltay = 0.50*(y(2)-y(1));
    V = A*deltay;
    
    for k=2:1:n_bins+1
        if (phi(i,1) < bins(k))
            sum_bins(k-1,1) = sum_bins(k-1,1)+V;
            break;
        end
    end
    
end

for j=2:1:ny-1
    
    hm = y(j) - 0.5*(y(j)-y(j-1));
    hp = y(j) + 0.5*(y(j+1)-y(j));
    deltay = hp-hm;
    
    %West
    ri = x(1);
    re = x(1)+0.50*(x(2)-x(1));
    A = pi*(re^2-ri^2);
    V = A*deltay;
    for k=2:1:n_bins+1
        if (phi(1,j) < bins(k))
            sum_bins(k-1,1) = sum_bins(k-1,1)+V;
            break;
        end
    end
    
    %East
    ri = x(nx-1)+0.50*(x(nx)-x(nx-1));
    re = x(nx);
    A = pi*(re^2-ri^2);
    V = A*deltay;
    for k=2:1:n_bins+1
        if (phi(nx,j) < bins(k))
            sum_bins(k-1,1) = sum_bins(k-1,1)+V;
            break;
        end
    end
    
end

%Corners S
deltay = 0.50*(y(2)-y(1));

%W
ri = x(1);
re = x(1)+0.50*(x(2)-x(1));
A = pi*(re^2-ri^2);
V = A*deltay;
for k=2:1:n_bins+1
    if (phi(1,1) < bins(k))
        sum_bins(k-1,1) = sum_bins(k-1,1)+V;
        break;
    end
end
%E
ri = x(nx-1)+0.50*(x(nx)-x(nx-1));
re = x(nx);
A = pi*(re^2-ri^2);
V = A*deltay;
for k=2:1:n_bins+1
    if (phi(nx,1) < bins(k))
        sum_bins(k-1,1) = sum_bins(k-1,1)+V;
        break;
    end
end

%Corners N
deltay = 0.50*(y(nx)-y(nx-1));

%W
ri = x(1);
re = x(1)+0.50*(x(2)-x(1));
A = pi*(re^2-ri^2);
V = A*deltay;
for k=2:1:n_bins+1
    if (phi(1,ny) < bins(k))
        sum_bins(k-1,1) = sum_bins(k-1,1)+V;
        break;
    end
end
%E
ri = x(nx-1)+0.50*(x(nx)-x(nx-1));
re = x(nx);
A = pi*(re^2-ri^2);
V = A*deltay;
for k=2:1:n_bins+1
    if (phi(nx,ny) < bins(k))
        sum_bins(k-1,1) = sum_bins(k-1,1)+V;
        break;
    end
end

sum = 0;
for k=1:1:n_bins-1
    sum = sum + sum_bins(k,1);
end

if (planar == true)
    Vtot = (x(nx)-x(1))*(y(ny)-y(1));
else
    Vtot = pi*(x(nx)^2-x(1)^2)*(y(ny)-y(1));
end

if ( abs(Vtot-sum) > 0.001*Vtot)
    error('Error in pore size distribution');
end

for k=1:1:n_bins-1
    bin_y(k) = sum_bins(k,1)/sum;
end

for k=2:1:n_bins
    bin_x(k-1) = 0.5*(bins(k-1)+bins(k));
end
