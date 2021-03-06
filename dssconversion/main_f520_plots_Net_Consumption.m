clear all
close all
clc

%exctracts data from the figures obtained from the daily simulation runs (DailySimulations.m)
%extracts total MW with and without pv. Also calculates total pv generation in MW from the .mat 
%files generated by the main_f520_scenario1.m, main_f520_scenario2.m,
%main_f520_scenario3.m, and main_f520_scenario4.m


strDirectory='C:\Work\Projects\2012\1787-UCSD_PV\Simulation\System\SDGE';
cd(strDirectory); 
curdir = pwd;
for i=1:4
    switch i
        case 1
            path = [curdir '/tmp_scnr1'];
            fn = dir([path '/scenario1_*.mat']);
            for j=0:6  
                result_dir = [path '/f520_case' num2str(j) '/results'];
                SimName = ['Case ' num2str(j)];
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                load([path '/' fn(j+1).name]);
                o = [];
                for i = 1:length(c.pvsystem)
                    cn = c.pvsystem(i).daily;
                    [id, id] = ismember(cn,{c.loadshape.Name}');
                    if i == 1
                        o = c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    else
                        o = o + c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    end
                end
                if length(o) ==2880
                    xrange = 1:2880;
                    xrange = xrange./120;
%                     dt = 30/3600/24;
                elseif length(o) == 24
                    xrange = 1:24;
%                     dt = 1/24;
                else
                    error('error');
                end
                f = figure; hold on;
                h = plot(xrange, MWPV, 'r', xrange, MWNoPV, 'g', xrange, o/1000, 'b');
                legend(h,{'MW (PV)', 'MW (No PV)', 'MW (PV Output)'},'fontsize',20);
                title('Net Consumption and PV Output','fontsize',20)
                xlabel('Time of Day','fontsize',20)
                ylabel('Power','fontsize',20)
                grid on;
                xlim([min(xrange) max(xrange)])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
                close all;
            end
        case 2
            disp('________________________________________________________________');
            path = [curdir '/tmp_scnr2'];
            fn = dir([path '/scenario2_*.mat']);
            for j=0:6
                result_dir = [path '/f520_case' num2str(j) '/results'];
                SimName = ['Case ' num2str(j)];
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                load([path '/' fn(j+1).name]);
                o = [];
                for i = 1:length(c.pvsystem)
                    cn = c.pvsystem(i).daily;
                    [id, id] = ismember(cn,{c.loadshape.Name}');
                    if i == 1
                        o = c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    else
                        o = o + c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    end
                end
                if length(o) ==2880
                    xrange = 1:2880;
                    xrange = xrange./120;
%                     dt = 30/3600/24;
                elseif length(o) == 24
                    xrange = 1:24;
%                     dt = 1/24;
                else
                    error('error');
                end
                f = figure; hold on;
                h = plot(xrange, MWPV, 'r', xrange, MWNoPV, 'g', xrange, o/1000, 'b');
                legend(h,{'MW (PV)', 'MW (No PV)', 'MW (PV Output)'},'fontsize',20);
                title('Net Consumption and PV Output','fontsize',20)
                xlabel('Time of Day','fontsize',20)
                ylabel('Power','fontsize',20)
                grid on;
                xlim([min(xrange) max(xrange)])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
                close all;
            end
        case 3
            disp('________________________________________________________________');
            path = [curdir '/tmp_scnr3'];
            fn = dir([path '/scenario3_*.mat']);
            for j=0:6
                result_dir = [path '/f520_case' num2str(j) '/results'];
                SimName = ['Case ' num2str(j)];
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                load([path '/' fn(j+1).name]);
                o = [];
                for i = 1:length(c.pvsystem)
                    cn = c.pvsystem(i).daily;
                    [id, id] = ismember(cn,{c.loadshape.Name}');
                    if i == 1
                        o = c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    else
                        o = o + c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    end
                end
                if length(o) ==2880
                    xrange = 1:2880;
                    xrange = xrange./120;
%                     dt = 30/3600/24;
                elseif length(o) == 24
                    xrange = 1:24;
%                     dt = 1/24;
                else
                    error('error');
                end
                f = figure; hold on;
                h = plot(xrange, MWPV, 'r', xrange, MWNoPV, 'g', xrange, o/1000, 'b');
                legend(h,{'MW (PV)', 'MW (No PV)', 'MW (PV Output)'},'fontsize',20);
                title('Net Consumption and PV Output','fontsize',20)
                xlabel('Time of Day','fontsize',20)
                ylabel('Power','fontsize',20)
                grid on;
                xlim([min(xrange) max(xrange)])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
                close all;
            end
       case 4
            disp('________________________________________________________________');
            path = [curdir '/tmp_scnr4'];
            fn = dir([path '/scenario4_*.mat']);
            for j=0:6
                result_dir = [path '/f520_case' num2str(j) '/results'];
                SimName = ['Case ' num2str(j)];
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                load([path '/' fn(j+1).name]);
                o = [];
                for i = 1:length(c.pvsystem)
                    cn = c.pvsystem(i).daily;
                    [id, id] = ismember(cn,{c.loadshape.Name}');
                    if i == 1
                        o = c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    else
                        o = o + c.loadshape(id).Mult*c.pvsystem(i).kVA;
                    end
                end
                if length(o) ==2880
                    xrange = 1:2880;
                    xrange = xrange./120;
%                     dt = 30/3600/24;
                elseif length(o) == 24
                    xrange = 1:24;
%                     dt = 1/24;
                else
                    error('error');
                end
                f = figure; hold on;
                h = plot(xrange, MWPV, 'r', xrange, MWNoPV, 'g', xrange, o/1000, 'b');
                legend(h,{'MW (PV)', 'MW (No PV)', 'MW (PV Output)'},'fontsize',20);
                title('Net Consumption and PV Output','fontsize',20)
                xlabel('Time of Day','fontsize',20)
                ylabel('Power','fontsize',20)
                grid on;
                xlim([min(xrange) max(xrange)])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
                saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
                close all;
            end
    end
end