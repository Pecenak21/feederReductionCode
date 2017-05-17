clear
clc

pathToFile='C:\Users\Zactus\FeederReductionRepo\8500-Node\Master.dss';
% pathToFile='C:\Users\Zactus\FeederReductionRepo\ckt5\Master_ckt5.dss';
% pathToFile='c:\users\zactus\FeederReductionRepo\123Bus\IEEE123Master.dss';

	o = actxserver('OpendssEngine.dss');
	dssText = o.Text; dssText.Command = 'Clear';
	dssText.Command = ['Compile "' pathToFile '"'];
	dssCircuit = o.ActiveCircuit;
	circuit.buslist.id=regexprep(dssCircuit.AllBUSNames,'-','_');
	buslist=circuit.buslist.id;
	circuit.buslist.coord=zeros(length(circuit.buslist.id),2);
	delete(o);
% 	criticalBuses=circuit.buslist.id;%([1,100,1000]);
	kill=find(~cellfun(@isempty,regexp(buslist,'x')));
	kill(1:4)=[];
	buslist(kill)=[]; criticalBuses=buslist;
	cd c:/users/zactus/FeederReductionRepo/
[circuit, circuit_orig] = reducingFeeders(pathToFile,criticalBuses)