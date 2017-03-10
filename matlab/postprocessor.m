clear all;

%folder name
planar_symmetry = false;
list_times = [6 12 18 24 30 36 48 54 60 66 72];
list_disks = [1];
folder_name = 'C:\Users\acuoci\Desktop\Brembo-Simulations\AnlageA\Densification\SensitivityAnalysesJanuary2017\Baseline\Sim31';

%Number of times
list_times_size = size(list_times);
n_times = list_times_size(2);

%Number of disks
list_disks_size = size(list_disks);
n_disks = list_disks_size(2);

%Loop
for k=1:1:n_disks

    %reconstruct disk file name
    disk_file_name = strcat(folder_name,'\Output.', int2str(list_disks(k)));
    
    fprintf('Processing disk %d ...\n', k);
    
    %Disk
    disk(k) = DiskClass(n_times, planar_symmetry);
    
    for time=1:1:n_times

        %reconstruct file name
        file_name = strcat(disk_file_name, '\Solution.', int2str(list_times(time)), '.out');
        
        %open file
        fprintf(' - Reading time %d ...\n', time);
        
        if (time == 1)
            disk(k) = disk(k).readFromFileFirstTime(file_name);
        end
        
        disk(k) = disk(k).readFromFile(file_name, list_times(time));
            
        disk(k) = disk(k).calculateMeanFields(time);
        disk(k) = disk(k).calculateMinFields(time);
        disk(k) = disk(k).calculateMaxFields(time);
        disk(k) = disk(k).calculateSampleProfiles(time);
        
    end
end

f = figure;
set(f,'name','Figure: Mean C2H2 mole fraction','numbertitle','off')
hold on;
for k=1:1:n_disks
    display_name = strcat('Disk ', int2str(k));
    plot (list_times, disk(k).mean_x_C2H2, '-o', 'LineWidth',2, 'DisplayName', display_name);
end
title('Mean C2H2 mole fraction');
xlabel('time [h]');
ylabel('mole fraction');
legend('show');
hold off;

f = figure;
set(f,'name','Figure: Mean C6H6 mole fraction','numbertitle','off')
hold on;
for k=1:1:n_disks
    display_name = strcat('Disk ', int2str(k));
    plot (list_times, disk(k).mean_x_C6H6, '-o', 'LineWidth',2, 'DisplayName', display_name);
end
title('Mean C6H6 mole fraction');
xlabel('time [h]');
ylabel('mole fraction');
legend('show');
hold off;

f = figure;
set(f,'name','Figure: Mean C10H8 mole fraction','numbertitle','off')
hold on;
for k=1:1:n_disks
    display_name = strcat('Disk ', int2str(k));
    plot (list_times, disk(k).mean_x_C10H8, '-o', 'LineWidth',2, 'DisplayName', display_name);
end
title('Mean C10H8 mole fraction');
xlabel('time [h]');
ylabel('mole fraction');
legend('show');
hold off;

f = figure;
set(f,'name','Figure: Mean H2 mole fraction','numbertitle','off')
hold on;
for k=1:1:n_disks
    display_name = strcat('Disk ', int2str(k));
    plot (list_times, disk(k).mean_x_H2, '-o', 'LineWidth',2, 'DisplayName', display_name);
end
title('Mean H2 mole fraction');
xlabel('time [h]');
ylabel('mole fraction');
legend('show');
hold off;

f = figure;
set(f,'name','Figure: Mean C2H2/C6H6 ratio', 'numbertitle','off')
hold on;
for k=1:1:n_disks
    display_name = strcat('Disk ', int2str(k));
    plot (list_times, disk(k).mean_x_C2H2./disk(k).mean_x_C6H6, '-o', 'LineWidth',2, 'DisplayName', display_name);
end
title('Mean C2H2/C6H6 ratio');
xlabel('time [h]');
ylabel('mole fraction');
legend('show');
hold off;

for time=1:1:n_times

    a = [];
    for k=1:1:n_disks
        a=[a disk(k).bin_y(time,:)'];
    end
    
    f = figure;
    set(f,'name','Pore size distribution', 'numbertitle','off')
    hold on;
    bar( disk(1).bin_x(time,:), a);
    plot_title = strcat('Pore size distribution @ ', int2str(list_times(time)), ' h'); 
    title(plot_title);
    xlabel('pore radius [micron]');
    ylabel('percentage [%]');
    legend('show');
    hold off;
    
end
    
