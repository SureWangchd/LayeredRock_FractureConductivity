function info=r7_ct_native_solver(m,tol)
% Retain COMSOL's automatic stationary sequence, but forbid a fine-grid
% direct factorization at ~2.2 million DOFs. Use native GMRES + native AMG.
if nargin<2,tol=1e-10;end
m.study('std1').feature('stat').set('usestol',true);
m.study('std1').feature('stat').set('stol',tol);
validateattributes(tol,{'numeric'},{'scalar','real','finite','positive','<',1});
tags=cell(m.sol.tags);if isempty(tags),m.study('std1').createAutoSequences('all');end
st=m.sol('sol1').feature('s1');st.set('control','user');st.set('stol',tol);
children=cell(st.feature.tags);iterTag='';fcTags={};
for j=1:numel(children)
 type=char(st.feature(children{j}).getType);
 if strcmp(type,'Iterative') && isempty(iterTag),iterTag=children{j};end
 if strcmp(type,'FullyCoupled'),fcTags{end+1}=children{j};end %#ok<AGROW>
end
if isempty(iterTag),iterTag='iCT';st.create(iterTag,'Iterative');end
it=st.feature(iterTag);it.set('linsolver','gmres');it.set('itrestart',50);it.set('maxlinit',2000);
sub=cell(it.feature.tags);mgTag='';
for j=1:numel(sub),if strcmp(char(it.feature(sub{j}).getType),'Multigrid'),mgTag=sub{j};break;end,end
if isempty(mgTag),mgTag='mgCT';it.create(mgTag,'Multigrid');end
mg=it.feature(mgTag);mg.set('prefun','amg');mg.set('amgauto',3);mg.set('mglevels',8);mg.set('maxcoarsedof',5000);
% An automatically created ILU or other sibling preconditioner must not
% silently take precedence over the explicitly selected native AMG node.
sub=cell(it.feature.tags);for j=1:numel(sub),it.feature(sub{j}).active(strcmp(sub{j},mgTag));end
if isempty(fcTags),st.create('fcCT','FullyCoupled');fcTags={'fcCT'};end
for j=1:numel(fcTags),st.feature(fcTags{j}).set('linsolver',iterTag);end
actual=struct();unavailable=struct();
props={'rhob','prefuntype','irestol','usenlweights'};
for j=1:numel(props)
 name=props{j};
 try
  if strcmp(name,'prefuntype'),value=char(it.getString(name));
  elseif strcmp(name,'usenlweights'),value=logical(it.getBoolean(name));
  else,value=double(it.getDouble(name));end
  actual.(name)=value;
 catch ME
  actual.(name)='unavailable';unavailable.(name)=ME.message;
 end
end
info=struct('implementation','Native COMSOL automatic stationary sequence with explicit GMRES+AMG selection; no external solver',...
 'relative_tolerance',tol,'linear_solver','gmres','restart',50,'max_linear_iterations',2000,...
 'preconditioner','COMSOL algebraic multigrid','AMG_quality',3,'max_levels',8,'max_coarse_dof',5000,...
 'stationary_feature','s1','iterative_feature',iterTag,'multigrid_feature',mgTag,...
 'fine_grid_direct_factorization_allowed',false,'fully_coupled_features',{fcTags},...
 'actual_iterative_properties',actual,'unavailable_property_details',unavailable,...
 'verification_note','API properties checked against local COMSOL6.3 docs; actual convergence and conservation must be verified for each state.');
end

