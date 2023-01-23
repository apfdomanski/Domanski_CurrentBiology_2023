function save_parfor(fname, data)
% Permits saving by functions whilst operating within a parfor loop
var_name=genvarname(inputname(2));
eval([var_name '=data;']);

try
    save(fname,var_name,'-append');
catch
    save(fname,var_name);
end