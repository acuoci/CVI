classdef DiskClass
   
   %Properties
   properties
   
        A
        nc
        npoints
        variables_names
        
        index_time
        index_temperature
        index_porosity
        index_bulk_density
        index_porous_radius
        index_x_coordinate
        index_y_coordinate
        index_x_C2H2
        index_x_C6H6
        index_x_C10H8
        index_x_H2
        
        n_times
        nx
        ny
        x
        y
        
        porosity
        bulk_density
        porous_radius
        x_C2H2
        x_C6H6
        x_C10H8
        x_H2
        
        mean_porosity
        mean_bulk_density
        mean_porous_radius
        mean_x_C2H2
        mean_x_C6H6
        mean_x_C10H8
        mean_x_H2

        max_bulk_density
        max_porosity
        max_porous_radius
        max_x_C2H2
        max_x_C6H6
        max_x_C10H8
        max_x_H2
        
        min_bulk_density
        min_porosity
        min_porous_radius
        min_x_C2H2
        min_x_C6H6
        min_x_C10H8
        min_x_H2
        
        sample_x_center_bulk_density
        sample_y_center_bulk_density
        
        bin_x
        bin_y
       
        planar_symmetry
        current_time
   end
   
   %Methods
   methods
       
       function obj = DiskClass(n_times, planar_symmetry)
           
           obj.planar_symmetry = planar_symmetry;
           obj.n_times = n_times;
           obj = obj.memoryPreallocationTimeVariables;
       
       end
      
       function obj = readFromFileFirstTime(obj, file_name)
           
            %Open file
            fid = fopen(file_name,'r');

            %Read the first line of the file containing the header information
            headerline  = fgetl(fid);
            obj.variables_names = strsplit(headerline);
            variables_names_size = size(obj.variables_names);
            obj.nc = variables_names_size(2)-1;

            %Count the number of points
            obj.npoints = 0;
            while ~feof(fid)
                obj.npoints = obj.npoints+1;
                line = fgetl(fid);
            end
            
            %Memory allocation
            obj.A = zeros(obj.npoints,obj.nc);
            
            %Close the file
            fclose(fid);
            
            %Read data
            obj = obj.readFromFile(file_name, 0);
            
            %Recognize indices
            obj = obj.recognizeIndices;
            
            %Recognize coordinates
            obj = obj.recognizeCoordinates;
            
            %Memory preallocation
            obj = obj.memoryPreallocation;
          
       end
       
       function obj = readFromFile(obj, file_name, current_time)
           
            %Current time
            obj.current_time = current_time;
            
            %Open file
            fid = fopen(file_name,'r');

            %Read data (header line)
            headerline = fgetl(fid);
            
            %Read data
            for j=1:1:obj.npoints
                line = fgetl(fid);
                line_elements = strsplit(line);
                for i=1:1:obj.nc
                    obj.A(j, i) = str2double(line_elements(1,i));
                end
            end
            
            %Close the file
            fclose(fid);
            
            %Fill fields
            obj = obj.fillFields;
            
       end
       
       function obj = recognizeIndices(obj)
       
            for i=1:1:obj.nc
                
                if (strncmpi('time[s]',obj.variables_names(i),7) == true)
                    obj.index_time = i;
                end
                if (strncmpi('T[K]',obj.variables_names(i),4) == true)
                    obj.index_temperature = i;
                end
                if (strncmpi('eps[-]',obj.variables_names(i),6) == true)
                    obj.index_porosity = i;
                end 
                if (strncmpi('rhoB',obj.variables_names(i),4) == true)
                    obj.index_bulk_density = i;
                end 
                if (strncmpi('rp[micron]',obj.variables_names(i),10) == true)
                    obj.index_porous_radius = i;
                end 
                if (strncmpi('x[mm]',obj.variables_names(i),5) == true)
                    obj.index_x_coordinate = i;
                end
                if (strncmpi('y[mm]',obj.variables_names(i),5) == true)
                    obj.index_y_coordinate = i;
                end
                if (strncmpi('C2H2_w',obj.variables_names(i),6) == true)
                    obj.index_x_C2H2 = i;
                end
                if (strncmpi('C6H6_w',obj.variables_names(i),6) == true)
                    obj.index_x_C6H6 = i;
                end
                if (strncmpi('C10H8_w',obj.variables_names(i),7) == true)
                    obj.index_x_C10H8 = i;
                end
                if (strncmpi('H2_w',obj.variables_names(i),4) == true)
                    obj.index_x_H2 = i;
                end

            end
       end
       
       function obj = recognizeCoordinates(obj)
       
            obj.nx = 1;
            obj.x(1) = obj.A(1, obj.index_x_coordinate);
            for i=2:1:obj.npoints
                if (obj.A(i, obj.index_x_coordinate) == obj.A(1, obj.index_x_coordinate))
                    break;
                else
                    obj.x(i) = obj.A(i, obj.index_x_coordinate);
                    obj.nx = obj.nx+1;
                end
            end

            obj.ny = obj.npoints/obj.nx;
            for i=1:1:obj.ny
                obj.y(i) = obj.A((i-1)*obj.nx+1, obj.index_y_coordinate);
            end
            
       end
       
       function obj = memoryPreallocation(obj)
           
            obj.porosity = zeros(obj.nx,obj.ny);
            obj.bulk_density = zeros(obj.nx,obj.ny);
            obj.porous_radius = zeros(obj.nx,obj.ny);
            obj.x_C2H2 = zeros(obj.nx,obj.ny);
            obj.x_C6H6 = zeros(obj.nx,obj.ny);
            obj.x_C10H8 = zeros(obj.nx,obj.ny);
            obj.x_H2 = zeros(obj.nx,obj.ny);
            
            obj.sample_x_center_bulk_density = zeros(obj.n_times, obj.nx);
            obj.sample_y_center_bulk_density = zeros(obj.n_times, obj.ny);
            
       end
       
       function obj = memoryPreallocationTimeVariables(obj)
           
            obj.mean_bulk_density = zeros(obj.n_times,1);
            obj.mean_porosity = zeros(obj.n_times,1);
            obj.mean_porous_radius = zeros(obj.n_times,1);
            obj.mean_x_C2H2 = zeros(obj.n_times,1);
            obj.mean_x_C6H6 = zeros(obj.n_times,1);
            obj.mean_x_C10H8 = zeros(obj.n_times,1);
            obj.mean_x_H2 = zeros(obj.n_times,1);
            
            obj.max_bulk_density = zeros(obj.n_times,1);
            obj.max_porosity = zeros(obj.n_times,1);
            obj.max_porous_radius = zeros(obj.n_times,1);
            obj.max_x_C2H2 = zeros(obj.n_times,1);
            obj.max_x_C6H6 = zeros(obj.n_times,1);
            obj.max_x_C10H8 = zeros(obj.n_times,1);
            obj.max_x_H2 = zeros(obj.n_times,1);
            
            obj.min_bulk_density = zeros(obj.n_times,1);
            obj.min_porosity = zeros(obj.n_times,1);
            obj.min_porous_radius = zeros(obj.n_times,1);
            obj.min_x_C2H2 = zeros(obj.n_times,1);
            obj.min_x_C6H6 = zeros(obj.n_times,1);
            obj.min_x_C10H8 = zeros(obj.n_times,1);
            obj.min_x_H2 = zeros(obj.n_times,1);
            
       end
       
       function obj = fillFields(obj)
           
            k = 1;
            for j=1:1:obj.ny
                for i=1:1:obj.nx
                    obj.porosity(i,j) = obj.A(k, obj.index_porosity);
                    obj.bulk_density(i,j) = obj.A(k, obj.index_bulk_density);
                    obj.porous_radius(i,j) = obj.A(k, obj.index_porous_radius);
                    obj.x_C2H2(i,j) = obj.A(k, obj.index_x_C2H2);
                    obj.x_C6H6(i,j) = obj.A(k, obj.index_x_C6H6);
                    obj.x_C10H8(i,j) = obj.A(k, obj.index_x_C10H8);
                    obj.x_H2(i,j) = obj.A(k, obj.index_x_H2);
                    k = k+1;
                end
            end
            
            f = figure;
            set(f,'name','Figure: C2H2/C6H6 ratio', 'numbertitle','off')
            hold on;
            surf(obj.x, obj.y, obj.x_C2H2./obj.x_C6H6);
            plot_title = strcat('C2H2/C6H6 ratio @', int2str(obj.current_time), ' h');
            title(plot_title);
            xlabel('radial coordinate [mm]');
            ylabel('axial coordinate [mm]');
            colorbar;
            hold off;
            
       end
       
       function obj = calculateMeanFields(obj, time)
           
            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.bulk_density);
            obj.mean_bulk_density(time) = mean;
            
            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.porosity);
            obj.mean_porosity(time) = mean;
            
            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.porous_radius);
            obj.mean_porous_radius(time) = mean;
            
            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.x_C2H2);
            obj.mean_x_C2H2(time) = mean;

            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.x_C6H6);
            obj.mean_x_C6H6(time) = mean;
            
            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.x_C10H8);
            obj.mean_x_C10H8(time) = mean;

            [~, mean] = integral2D(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.x_H2);
            obj.mean_x_H2(time) = mean;
            
            [obj.bin_x(time,:), obj.bin_y(time,:)] = porousDistribution(obj.planar_symmetry, obj.nx, obj.x, obj.ny,obj.y, obj.porous_radius);
       
       end
       
       function obj = calculateMaxFields(obj, time)
           
            obj.max_bulk_density(time) = max(max(obj.bulk_density));
            
            obj.max_porosity(time) = max(max(obj.porosity));
            
            obj.max_porous_radius(time) = max(max(obj.porous_radius));
          
            obj.max_x_C2H2(time) = max(max(obj.x_C2H2));
            
            obj.max_x_C6H6(time) = max(max(obj.x_C6H6));
            
            obj.max_x_C10H8(time) = max(max(obj.x_C10H8));
            
            obj.max_x_H2(time) = max(max(obj.x_H2));
      
       end
       
       function obj = calculateMinFields(obj, time)

            obj.min_bulk_density(time) = min(min(obj.bulk_density));
            
            obj.min_porosity(time) = min(min(obj.porosity));
            
            obj.min_porous_radius(time) = min(min(obj.porous_radius));

            obj.min_x_C2H2(time) = min(min(obj.x_C2H2));

            obj.min_x_C6H6(time) = min(min(obj.x_C6H6));
            
            obj.min_x_C10H8(time) = min(min(obj.x_C10H8));
            
            obj.min_x_H2(time) = min(min(obj.x_H2));
       
       end
       
       function obj = calculateSampleProfiles(obj, time)
           
            j_center = round(obj.ny/2);
            obj.sample_x_center_bulk_density(time,:) = obj.bulk_density(:,j_center)';

            i_center = round(obj.nx/2);
            obj.sample_y_center_bulk_density(time,:) = obj.bulk_density(i_center,:);
        
       end
       
   end
   
end