
%%testvalue = "ODEs:\n     d(Rii)/dt = 1/spine*(-reaction_1 + reaction_2)\n     d(cAMP)/dt = 1/spine*(-reaction_5 - reaction_6)\n     d(RiiP)/dt = 1/spine*(-reaction_2 - reaction_4 - reaction_6)\n     d(Rii_C)/dt = 1/spine*(reaction_1 - reaction_3)\n     d(RiiP_cAMP)/dt = 1/spine*(reaction_6 - reaction_7)\n     d(RiiP_C)/dt = 1/spine*(reaction_3 + reaction_4 - reaction_5)\n     d(RiiP_C_cAMP)/dt = 1/spine*(reaction_5 + reaction_7)\n     d(C)/dt = 1/spine*(-reaction_1 - reaction_4 - reaction_7)\n     d(temp)/dt = rateout(time)\n     \n     Fluxes:\n     reaction_1 = (kf_RiiP_C*C*Rii)*spine-(kb_RiiP_C*Rii_C)*spine\n     reaction_2 = (kf_RiiP*RiiP)*spine\n     reaction_3 = (kf_Rii_C__RiiP_C*Rii_C)*spine-(kb_Rii_C__RiiP_C*RiiP_C)*spine\n     reaction_4 = (kf_RiiPxC__RiiP_C*RiiP*C)*spine-(kb_RiiPxC__RiiP_C*RiiP_C)*spine\n     reaction_5 = (kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP)*spine-(kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP)*spine\n     reaction_6 = (kf_cAMPxRiiP__Rii_cAMP*cAMP*RiiP)*spine-(kb_cAMPxRiiP__Rii_cAMP*RiiP_cAMP)*spine\n     reaction_7 = (kf_Rii_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C)*spine-(kb_Rii_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP)*spine\n     \n     Parameter Values:\n     kf_RiiP_C = 0.52\n     kb_RiiP_C = 0.0003\n     kf_RiiP = 0.032\n     kf_Rii_C__RiiP_C = 0.022\n     kf_RiiP_CxcAMP__RiiP_C_cAMP = 0.496\n     kb_RiiP_CxcAMP__RiiP_C_cAMP = 1.413\n     kf_Rii_cAMPxC__RiiP_C_cAMP = 0.07\n     kb_Rii_cAMPxC__RiiP_C_cAMP = 0.1046\n     kb_cAMPxRiiP__Rii_cAMP = 1.413\n     kb_Rii_C__RiiP_C = 0\n     kf_cAMPxRiiP__Rii_cAMP = 5.5\n     kf_RiiPxC__RiiP_C = 1.26\n     kb_RiiPxC__RiiP_C = 0.018\n     spine = 1\n     \n     Initial Conditions:\n     Rii = 13.704\n     cAMP = 5\n     RiiP = 0.297\n     Rii_C = 0.437\n     RiiP_cAMP = 0\n     RiiP_C = 0.5616\n     RiiP_C_cAMP = 0\n     C = 0.001365\n     temp = -1";

function [F,M]=make_flux_function(SimbiologyEquations,varargin)
  %%
  %% Usage: [F,M]=make_flux(SimbiologyEquations,[file_name])
  %%
  %% SimbiologyEquations is the text that is returned by simbiology's
  %% "getequations".
  %% If no output arguments are specified, then the
  %% flux function is written into a file (inside the current working
  %% directory).
  %% The filename will be automatically generated; an
  %% optional second argument to this function (file_name) will
  %% override this.
  
  %% transform the one character array into a cell array of strings
  eq=strsplit(SimbiologyEquations,'\n');
  %% remove all empty lines
  l=cellfun(@isempty,regexp(eq,'^\s*$','start','emptymatch'));
  eq=regexprep(eq(l),'\s',''); % all spaces and tabs removed
  %% find all Headers:
  l=contains(eq,':');
  %% alternatively:
  %%s=strfind(eq,":");
  %%l=~cellfun(@isempty,s);
  H=regexprep(eq(l),':','');% these are all the headers, like "ODEs"
  nH=length(H);
  assert(nH>0);
  V=cell(1,nH);
  i=find(l); % these are the line-numbers of the Headers
  assert(length(i)>0);
  N=diff([i,length(eq)]); % numbers of lines between headers
  for j=1:length(H)
    n=N(j)-1;
    V{j}={eq(i(j)+[1:n])};
  end%for
  HV=cat(1,H,V);
  M=struct(HV{:});
  %% make a function handle that maps a parameter vector p to fluxes:\
  if isfield(M,'ParameterValues') && isfield(M,'Fluxes') && isfield(M,'ODEs')
    np=length(M.ParameterValues);
    %% matlab has no «ostrsplit» function
    p=strsplit(sprintf('p(%i)\t',[1:np]),'\t');
    p=p(1:np);
    nx=length(M.ODEs);
    x=strsplit(sprintf('x(%i)\t',[1:nx]),'\t');
    x=x(1:nx);
    M.ParameterNames=regexprep(M.ParameterValues,'=.*$','');
    XNames=regexprep(M.ODEs,'^d\((\w+)\)/dt.*$','$1');
    M.StateVariableNames=XNames;                      % return just the match ... of "d(...)/dt"
    PPattern=strcat('\<',M.ParameterNames,'\>');     % so we get exact matches
    XPattern=strcat('\<',M.StateVariableNames,'\>'); % same
    Flux=regexprep(M.Fluxes,PPattern,p);      % replace parameter names with p(i)
    Flux=regexprep(Flux,XPattern,x);      % replace species names with x(i)
    M.Flux=regexprep(Flux,'^\w+=','');    
  else
    M
    error('ODEs, ParameterValues and Fluxes must be contained in the equations string');
  end%if
  F=@(time,x,p) [cellfun(@eval,M.Flux)]';
  if (nargout==0)
    if (nargin<2)
      fname=sprintf('flux_of_%s_model_with_%i_parameters',M.StateVariableNames{1},np);
    else
      fname=varargin{2};
    end%if
    fid=fopen(strcat(fname,'.m'),'w');
    fprintf(fid,'function [F]=%s(time,x,p)\n',fname);
    fprintf(fid,' %%%% Usage [F]=%s(time,x,p)\n',fname);
    fprintf(fid,' %%%% x: %i state variable values (species)\n',nx);
    fprintf(fid,' %%%% p: %i parameter values (e.g. kinetic coefficients)\n',np);
    fprintf(fid,' %%%% state variables:\n');
    fprintf(fid,' %%%%  %s\n',M.StateVariableNames{:});
    fprintf(fid,' %%%% parameters:\n');
    fprintf(fid,' %%%%  %s\n',M.ParameterNames{:});
    fprintf(fid,'\n');
    fprintf(fid,' assert(length(x)==%i);\n',nx);
    fprintf(fid,' assert(length(p)==%i);\n',np);
    fprintf(fid,' F=[...\n');
    fprintf(fid,'    %s;...\n',M.Flux{:});
    fprintf(fid,'   ];\n');
    fprintf(fid,'end%%function');
    fclose(fid);
  end%if
end%function

