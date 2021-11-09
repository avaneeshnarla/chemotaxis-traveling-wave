%% Initialise
all_files_info=dir('**/parameters.txt');
file_number=1;
data_table=table;
locn={};
%% Fill Table
for file_number=1:length(all_files_info)
    % Gather Expansion Speed
    max_loc_file=dir([all_files_info(file_number).folder,'/max_vel_loc*']);
    if length(max_loc_file)~=1
        continue;
    end
    load([max_loc_file.folder,'/',max_loc_file.name]);
    peak_loc_file=dir([all_files_info(file_number).folder,'/peak_loc*']);
    load([peak_loc_file.folder,'/',peak_loc_file.name]);
    x_maxvelloc=x(max_vel_loc>0);
    max_vel_loc=max_vel_loc(max_vel_loc>0);
    [max_height,max_time]=max(max_vel_loc(1001:end));
    peak_indx=max(max_time)+1000;
    if max(x_maxvelloc)>2.5e4
        x_maxvelloc=x_maxvelloc/60.0;
    end
    window_max=length(x_maxvelloc)/23;
    x1=x_maxvelloc(peak_indx-100-window_max:peak_indx-100);
    max_vel_loc1=max_vel_loc(peak_indx-100-window_max:peak_indx-100);
    [p,S]=polyfit(x1(max_vel_loc1>1),max_vel_loc1(max_vel_loc1>1),1);
    R2=(S.normr/norm(max_vel_loc1 - mean(max_vel_loc1)))^2;
    display(R2)
    if (R2>2e-3 || S.normr==0)
       continue
    end
    
    
    % Gather peak max and trough
    max_time=0;
    profile_num=0;
    list=dir([all_files_info(file_number).folder,'/density_profile*']);
    list2=dir([all_files_info(file_number).folder,'/attractant_profile*']);
    list3=dir([all_files_info(file_number).folder,'/velocity*']);
    
    for list_num=1:length(list)
        list_split=split(list(list_num).name,'_');
        last_num=list_split{length(list_split)};
        time_print=split(last_num,'.');
        if str2double(time_print{1})>max_time+100
            profile_num=list_num;
            max_time=str2double(time_print{1});
            dt=mod(str2double(time_print{1}),100);
        end
    end
    if profile_num==0
        continue
    end
    
    max_rho=0;
    min_rho=0;
    zmax=1;
    zmin=1;
    a_max=0;
    astar=0;
    load([list(profile_num).folder,'/',list(profile_num).name])
    load([list(profile_num).folder,'/',list2(profile_num).name])
    load([list(profile_num).folder,'/',list3(profile_num).name])
%     if ~length(findpeaks(rho_profile(1:3:end)))==1
%         continue
%     end
%     if ~isempty(findpeaks(a_profile(1:3:end),'MinPeakProminence',1e-3))
%         continue
%     end
    for index=4:length(rho_profile)-3
        if rho_profile(index)>rho_profile(index+3) && rho_profile(index)>rho_profile(index-3) && rho_profile(index)>max_rho
            max_rho=rho_profile(index);
            a_max=a_profile(index);
            zmax=index;
        elseif rho_profile(index)<rho_profile(index+3) && rho_profile(index)<rho_profile(index-3) && rho_profile(index)>1e-2
            min_rho=rho_profile(index);
            astar=a_profile(index);
            zmin=index;
        end
    end
    
    % Feed the table
    par_text=fopen([all_files_info(file_number).folder,'/parameters.txt']);
    txt=textscan(par_text,'%s');
    txt=txt{1,1};
    varTypes    = cell(1, 26);
    varTypes(1:26) = {'double'};
    varNames = {'T','dt','lmax','am','Da','D','a0','rho0','chi0','rhoc','ak'...
        ,'mu','ahill','gridsize','resln','exp_speed','max_rho','min_rho',...
        'min_a','zmax_zm','zmin_zm','amax','vmax','apmax','Ndiff','N0'};
    table_add = table('Size',[1 26],'VariableTypes',varTypes,'VariableNames',varNames);
    table_add.exp_speed=p(1);
    table_add.max_rho=max_rho;
    table_add.min_rho=min_rho;
    table_add.min_a=astar;
    locn{end+1}=all_files_info(file_number).folder;
    
    for txt_rows=1:length(txt)
        switch char(txt(txt_rows))
            case 'T'
                table_add.T=str2double(txt(txt_rows+2));
            case 'dt'
                %table_add.dt=str2double(txt(txt_rows+2));
                table_add.dt=dt;
            case 'lmax'
                table_add.lmax=str2double(txt(txt_rows+2));
            case 'am'
                if isnan(str2double(txt(txt_rows+2)))
                    dummy=char(txt(txt_rows+1));
                    table_add.am=str2double(dummy(2:end));
                else
                    table_add.am=str2double(txt(txt_rows+2));
                end
            case 'Da'
                table_add.Da=str2double(txt(txt_rows+2));
            case 'D'
                table_add.D=str2double(txt(txt_rows+2));
            case 'a_0'
                table_add.a0=str2double(txt(txt_rows+2));
            case 'rho0'
                table_add.rho0=str2double(txt(txt_rows+2));
            case 'chi_0'
                table_add.chi0=str2double(txt(txt_rows+2));
            case 'rhoc'
                table_add.rhoc=str2double(txt(txt_rows+2));
            case 'Kma'
                table_add.ak=str2double(txt(txt_rows+2));
            case 'ak'
                table_add.ak=str2double(txt(txt_rows+2));
            case 'alpha'
                table_add.mu=str2double(txt(txt_rows+2));
            case 'grid_size'
                table_add.gridsize=str2double(txt(txt_rows+2));
            case 'resln'
                table_add.resln=str2double(txt(txt_rows+2));
            case 'ahill'
                table_add.ahill=1;
        end
        row_text=char(txt(txt_rows));
        if length(row_text)>4 && all(row_text(1:4)=='rhoc')
            table_add.rhoc=str2double(row_text(5:end));
        elseif length(row_text)>3 && all(row_text(1:3)=='Kma')
            table_add.ak=str2double(row_text(4:end));
        elseif length(row_text)>5 && all(row_text(1:5)=='alpha')
            table_add.mu=str2double(row_text(6:end));
        elseif length(row_text)>5 && all(row_text(1:5)=='ahill')
            table_add.ahill=str2double(row_text(6:end));
        elseif length(row_text)>5 && all(row_text(1:5)=='chi_0')
            table_add.chi0=str2double(row_text(6:end));
        end
    end
    if table_add.gridsize==0
        table_add.gridsize=30;
    end
    [am_num,zm]=min(abs(a_profile-table_add.am));
    table_add.zmax_zm=(zmax-zm)*table_add.gridsize/length(a_profile);
    table_add.zmin_zm=(zmin-zm)*table_add.gridsize/length(a_profile);
    aprime=gradient(a_profile(1:3:end),table_add.gridsize/length(a_profile(1:3:end)));
    table_add.apmax=aprime(ceil(zmax/3));
    Ndi=trapz(rho_profile(zmax:3:end))*table_add.gridsize/length(rho_profile(1:3:end));
    table_add.Ndiff=Ndi;
    Ndi2=trapz(rho_profile(zmin:3:end))*table_add.gridsize/length(rho_profile(1:3:end));
    table_add.N0=Ndi2;
    table_add.amax=a_max;
    table_add.vmax=vel_profile(ceil(zmax/3));
    if iscell(vel_profile(ceil(zmax/3)))
        table_add.vmax=table_add.vmax{1};
    end
    fclose(par_text);
    %data_table=table('Size',[length(afi) 15],'VariableNames',{'T','dt','lmax','am','Da','D','a_0','rho0','chi_0','rhoc','Kma','alpha','ahill','grid_size','resln','exp_speed'});
    data_table=[data_table;table_add];
    %     if table_add.D==1e-3
    %     figure(99)
    %     clf;
    %     plot(x_maxvelloc,max_vel_loc,'o')
    %     hold on
    %     plot(x1,p(1)*x1+p(2),'linewidth',5)
    %     title({['chi=',num2str(table_add.chi0*1e6),' Da=',num2str(table_add.Da*1e6),' lmax=',num2str(table_add.lmax*3600),' D=',num2str(table_add.D*1e6),' Residual=',num2str(R2)],all_files_info(file_number).folder})
    %     saveas(gca,['exp_fits_1000/chi=',num2str(table_add.chi0*1e6),' Da=',num2str(table_add.Da*1e6),' lmax=',num2str(table_add.lmax*3600),' D=',num2str(table_add.D*1e6),'_',num2str(randi(100)),'.png'])
    %     end
end
[data_table1,unique_locn]=unique(data_table);
locn1=locn(1,unique_locn);
save('G:\My Drive\Projects\Hwa\Chemotaxis\Simulation Results\all_plots\all_data.mat','data_table1','locn1')