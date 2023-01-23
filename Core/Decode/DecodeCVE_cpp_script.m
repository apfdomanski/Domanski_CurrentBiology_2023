% DECODECVE_CPP_SCRIPT   Generate MEX-function DecodeCVE_mex from DecodeCVE.
% 
% Script generated from project 'DecodeCVE_cpp.prj' on 09-May-2018.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'DecodeCVE'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0);

%% Invoke MATLAB Coder.
cd('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\AssemblyCode\Core\Decode');
codegen -config cfg DecodeCVE -args ARGS{1} -nargout 1

