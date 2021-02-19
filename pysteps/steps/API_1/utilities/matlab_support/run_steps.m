
function [t, x, names] = run_steps(model, varargin)
    % Runs a Matlab Simbiology model using STEPS.
    %
    % model   a Matlab Simbiology model object.

    % varargin
    % - solver - 'wmdirect' or 'wmrk4' (deterministic), default: 'wmdirect' (stochastic),
    % - seed - set an explicit see for the pseudo random number generator. 
    %          If '0' or nothing given, a seed will be generated,
    % - realisations - the number of realisations to be computed, default: '1',
    % - rk4dt - internal time step of the Runge-Kutta solver wmrk4, default: '1e-5',
    % - dt - increment between the time points at which the molecular populations are stored,
    % - stoptime - the final time up to which the simulation will be computed,
    % - modeltofile - boolean flag, 'True' or False'. If set to 'True', 
    %                 the generated STEPS model is written to a file.

    % handle the varargs using Matlab's inputparser
    p = inputParser;

    % add name value pairs to the parsing scheme
    addParameter(p, 'solver', 'wmdirect', @(s) ismember(s, {'wmdirect', 'wmrk4'}));
    addParameter(p, 'seed', 0);
    addParameter(p, 'realisations', 1);
    addParameter(p, 'rk4dt', 1e-5);
    addParameter(p, 'dt', 0.1);
    addParameter(p, 'stoptime', 1.0);
    isb = @(f) isboolean(f);
    addParameter(p, 'modeltofile', false, @(f) islogical(f));
    
    parse(p, varargin{:})

    % show the parsed results
    disp(p.Results)

    % generate seed if required
    if(p.Results.seed == 0)
        rng_seed = uint32(randi([0, intmax('uint32')]));
    else
        rng_seed = p.Results.seed; 
    end

    % extract the model from Matlab by introspection
    [model_string, initial_condition_string] = extractModel(model);

    % get any events from the model
    event_str = getEvents(model);

    % create the reaction rate - reaction relation
    rate_names = char(rateToReaction(model));

    % call the Python side
    tofile = '';
    if(p.Results.modeltofile)
        tofile = {' --file '};
    end  
    
    % get the path of the matlab_support module
    steps_mod = py.importlib.import_module('steps.utilities.matlab_support');
    mod_path = char(steps_mod.abs_path);

    code = strcat({'python '}, mod_path, '/simulate.py ''', model_string, ''' ''', initial_condition_string, {''' '}, string(p.Results.stoptime), {' --events '}, {' ['''},event_str, ' , ' , rate_names, {'''] '},{' --seed '}, string(rng_seed), {' --solver '},  string(p.Results.solver), {' --dt '}, string(p.Results.dt), {' --rk4dt '}, string(p.Results.rk4dt),  {' --realisations '}, string(p.Results.realisations), tofile);

    [status, data] = system(char(code), '-echo');
    disp('STEPS exit status:')
    disp(status)

    % unwrap and return the json encoded data
    if(status ~= 0)
        print('Error occurred running STEPS: ' + data)
    end

    split_data = splitlines(data);
    out_data = split_data(end);
    result = jsondecode(char(out_data(1)));

    % generate the time vector
    % steps stores the molecular populations at
    % equidsitant points in time,
    % so there is no need to transfer the time vector
    t = linspace(0, p.Results.stoptime, size(result.data, 3));

    % check nargout to see
    % if we should return in the [t, x, names] convention
    % or as a SimData object 
    % we deviate from Matlab by allowing the [t, x, names] 
    % for multiple realisationsconvention
    if(nargout == 3)
        % t and names only once?
        x = cell.empty;
        names = result.species;
        disp('t, x, names output convention.')
       
        % iterate over all of the realisations
        for r=1:p.Results.realisations
            x{r} = [];
            names = [names; result.species];

            % iterate all the compartments
            for c=1:size(result.data, 2)
                x{r} = [x{r} squeeze(result.data(r,c,:,:))];
            end
        end

    elseif(nargout == 1)
        error('SimData object as output is not yet supported.')

        % ------------------------------------------------------
        % This branch is broken code and needs to be fixed.
        % It is rendered inactive by the previous error message.
        % ------------------------------------------------------
        
        % generate meta data Matlab is expecting
        simstr = struct('ModelName', model.name, 'ModelUUID', model.UUID);
        runstr = struct( ...
            'Variant', 'none', ...
            'ConfigSet', 'to be done', ...
            'SimulationDate', datetime, ...
            'SimulationType', 'STEPS'...
        );

        simDataObj = SimData.empty;
        for r=1:p.Results.realisations
            tmp = [];
            tmp_names = [];
            for c=1:size(result.data, 2)
                tmp = [tmp, squeeze(result.data(r,c,:,:))];
                tmp_names = [tmp_names; result.species];
            end
            % tmp = mat2dataset(squeeze(result.data(r,c,:,:)), 'VarNames', result.species);
            numel(num2cell(tmp,2))
            simDataObj(r) = SimData(t', tmp, simstr, runstr, num2cell(tmp,1))
        end
        % ----------------------
        % end of inactive code
        % ----------------------
    else
        error('Wrong number of output arguments: ', nargout, ' Expected either 1 or 3.')
    end
end


% recursive function getting all of the nested compartments
% try returning the array, no persistent variables

% get the compartment definition as well as the initial amounts
% in order to scale the reaction rate of zero order reactions 
% correctly, we also need to return the compartment's capacity
function [compartment_code, species_code, amount_code, capacity] = getCompartments(compartment, model, first_call)
    % persistent variable to hold the compartment data
    % we just append to this one
    % name, capacity
    list_of_compartments = zeros(2, 1); % not yet persistent

    % compartment
    name = compartment.get('Name');
    capacity = compartment.get('Capacity');
    unit = compartment.get('CapacityUnits');
        
    % check if the capacity unit is supported
    % NOTE: STEPS uses cubic meter for compartment volume
    factor = 1.0;
    if(strcmp(unit, 'l') || strcmp(unit, 'liter'))
        factor = 1e-3;
    elseif(strcmp(unit, 'ml') || strcmp(unit, 'milliliter'))
        factor = 1e-6;
    elseif(strcmp(unit, 'mul') || strcmp(unit, 'microliter'))
        factor = 1e-9;
    elseif(strcmp(unit, 'nl') || strcmp(unit, 'nanoliter'))
        factor = 1e-12;
    elseif(strcmp(unit, ''))
        warning('No capacity unit given, using default of 1 mul.')
        factor = 1e-9;
    else
        error('Unsupported capacity unit for compartment.')
    end    
    capacity = capacity * factor;
    
    % create the compartment in STEPS
    compartment_code = strcat(name, ' = swm.Comp("', name, '", geometry)\n');
    compartment_code = strcat(compartment_code, strcat(name, '.addVolsys("vsys")\n'));
    compartment_code = strcat(compartment_code, strcat(name, '.setVol(', num2str(capacity) ,')\n'));
    compartment_code = strcat(compartment_code, strcat('comp_names.append("', name, '")\n'));

    % extract and create the species in STEPS
    species_code = '';
    amount_code = '';
    species = compartment.get('Species');
    
    species_names = {};
    for s = species'
        if first_call & (isempty(find(ismember(species_names, s.get('Name')))))
            species_code = strcat(species_code, strcat(s.get('Name'), ' = smod.Spec("' , s.get('Name'), '", model)\n') );
            species_code = strcat(species_code, strcat('spec_names.append("', s.get('Name'), '")\n'));
            species_names = [species_names, s.get('Name')];
        end

        % we assume that Matlab takes the concentration as mol/l while STEPS assumes m^3
        if(strcmp(s.get('InitialAmountUnits'), 'molecule'))
            amount_code = strcat(amount_code, strcat('sim.setCompCount("', name, '", "', s.get('Name'), '", ', num2str(s.get('InitialAmount')), ')\n' ));
            amount_code = strcat(amount_code, strcat('comp_si_factor["', name, '"] = 1\n'));
        elseif(strcmp(s.get('initialAmountUnits'), 'micromolarity'))
            amount_code = strcat(amount_code, strcat('sim.setCompAmount("', name, '", "', s.get('Name'), '", ', num2str(s.get('InitialAmount') * 1e-3 * capacity), ')\n' ));
            amount_code = strcat(amount_code, strcat('comp_si_factor["', name, '"] = ', num2str(1e-3 * capacity), '\n'));
        elseif(strcmp(s.get('initialAmountUnits'), 'nanomolarity'))
            amount_code = strcat(amount_code, strcat('sim.setCompAmount("', name, '", "', s.get('Name'), '", ', num2str(s.get('InitialAmount') * 1e-6 * capacity), ')\n' ));
            amount_code = strcat(amount_code, strcat('comp_si_factor["', name, '"] = ', num2str(1e-6 * capacity), '\n'));
        elseif(strcmp(s.get('InitialAmountUnits'), ''))
            warning('No amount unit given, using default of molecules')
            amount_code = strcat(amount_code, strcat('sim.setCompCount("', name, '", "', s.get('Name'), '", ', num2str(s.get('InitialAmount')), ')\n' ));
            amount_code = strcat(amount_code, strcat('comp_si_factor["', name, '"] = 1\n'));
        else
            % raise exception, unkown unit
            error('Unknown unit.')
        end

        if(s.get('ConstantAmount'))
            amount_code = strcat(amount_code, strcat('sim.setCompClamped("', name, '", "', s.get('Name'), '", True)\n'));
        end
    end
    first_call = false;

    nested = compartment.get('Compartments');
    if(~isempty(nested))
        for n = nested
            tmp = getCompartments(n, model, first_call);
            cat(2, list_of_compartments, tmp);
        end
    end
end
% done getCompartments

function compartments = getReactionCompartments(r, model)
    compartments = {};
    tmp = strsplit(r.get('Reaction'), '->');
    for t = tmp
        tmp2 = strsplit(t{1}, '.');
        if(length(tmp2) == 2)
            compartments = [compartments, strtrim(tmp2(1))];
        else
            tmp = model.get('Compartments')
            if(length(tmp) ~= 1)
                error('Reaction definition does not contain compartment, so only one compartment in the model expected');
            end
            compartments = {tmp.get('Name')}; 
        end
    end
    
    % sort to preserve the order - compartments
    % a cell array
    compartments = sort(unique(compartments));
end

% handling functions
% patches is a cell array of patch names
function [code, patches] = handlePatch(compartments, patches)
    if(length(compartments) ~= 2)
        error('Expecting exactly 2 compartments.')
    end

    % arbitrary assignment
    inner_comp = compartments(1);
    outer_comp = compartments(2);

    % naming convention for patches: name of inner compartment _ name of outer compartment
    patch_name = strcat(inner_comp, '_', outer_comp);
 
    % check if patch already exists
    if(~isempty(strfind(patches, patch_name)))
        code = '';
        return
    end

    code = strcat(patch_name, ' = swm.Patch("', patch_name, '", geometry, ',  inner_comp, {', '}, outer_comp, ')\n');
    code = strcat(code, strcat(patch_name, '.addSurfsys("ssys")\n'));
    patches = unique([patches, patch_name]);
end

function [reaction_code, activation_code, sreac_to_patch] = handleSurfaceReaction(r, index, inner_comp, outer_comp, model, patches)
    reaction_code = '';
    activation_code = '';
    sreac_to_patch = '';

    eqn = r.get('Reaction');
    stoch = r.get('Stoichiometry');

    % parse the string - only consider non reversible reactions
    if(strfind(eqn, '<->'))
        error('Reversible reactions are not yet supported.')
    end

    tmp_code = strcat('r', num2str(index), ' = smod.SReac("r', num2str(index), '", ssys, ');
    tmp_olhs = 'olhs = [';
    tmp_orhs = 'orhs = [';
    tmp_ilhs = 'ilhs = [';
    tmp_irhs = 'irhs = [';

    i = 1;
    first = true;
    tmp = strsplit(r.get('Reaction'), '->');
    for t = tmp(1) % iterate reactants
        tmp2 = strtrim(strsplit(t{1}, '.'));
        sname = tmp2{2};
        comp = tmp2{1};

        for j = 1:abs(stoch(i))
            if(comp == inner_comp)
                if(first)
                    tmp_ilhs = strcat(tmp_ilhs, sname);
                    first = false;
                else
                    tmp_ilhs = strcat(tmp_ilhs, ', ', sname);
                end
            elseif(comp == outer_comp)
                if(first)
                    tmp_olhs = strcat(tmp_olhs, sname);
                    first = false;
                else
                    tmp_olhs = strcat(tmp_olhs, ', ', sname);
                end
            else
                error('Unknown compartment.')
            end
        end

        i = i + 1;
    end

    tmp_ilhs = strcat(tmp_ilhs, '], ');
    tmp_olhs = strcat(tmp_olhs, '], ');

    first = true;
    for t = tmp(2) % iterate products
        tmp2 = strtrim(strsplit(t{1}, '.'));
        comp = tmp2{1};
        sname = tmp2{2};

        for j = 1:abs(stoch(i))
            if(comp == inner_comp)
                if(first)
                    tmp_irhs = strcat(tmp_irhs, sname);
                    first = false;
                else
                    tmp_irhs = strcat(tmp_irhs, ', ', sname);
                end
            elseif(comp == outer_comp)
                if(first)
                    tmp_orhs = strcat(tmp_orhs, sname);
                    first = false;
                else
                    tmp_orhs = strcat(tmp_orhs, ', ', sname);
                end
            else
                error('Unknown compartment.')
            end
        end

        i = i + 1;
    end
    tmp_irhs = strcat(tmp_irhs, '], ');
    tmp_orhs = strcat(tmp_orhs, '], ');

    [rate_value, scaling_factor] = getRateValue(r, model, true);
    tmp_code = strcat(tmp_code, {' '}, tmp_ilhs, {' '}, tmp_irhs, {' '}, tmp_olhs, {' '}, tmp_orhs, ' kcst = ', num2str(rate_value), ' )\n' );        
    tmp_code = strcat(tmp_code, 'kcst_si_factor["r', num2str(index), '"] = ', num2str(scaling_factor), '\n');     

    reaction_code = strcat(reaction_code, strcat('rate', num2str(index), ' = ', num2str(getRateValue(r, model, true)), '\n'));
    reaction_code = strcat(reaction_code, tmp_code);

    reaction_code = strcat(reaction_code, strcat('sreac_to_patch["r', num2str(index) , '"] = "', inner_comp, '_', outer_comp, '"\n'));

    % deactivate in all patches
    for p = unique(patches)
        activation_code = strcat(activation_code, strcat('sim.setPatchSReacActive("', p, '", "r', num2str(index), '", False)\n'));
    end

    % activate only in the desired one
    activation_code = strcat(activation_code, strcat('sim.setPatchSReacActive("', inner_comp, '_', outer_comp, '", "r', num2str(index), '", True)\n'));
    
    activation_code = char(activation_code);
end

% build a json encoded string of all the reaction - reaction rate names
% this is an ordered list
function name = getRateName(reaction)
    if(isempty(reaction.get('KineticLaw')))
        if(isempty(reaction.get('ReactionRate')))
            error('Reaction rate missing: ', reaction)
        end
        name = reaction.get('ReactionRate');
    else
        name = reaction.get('KineticLaw').get('ParameterVariableNames');
    end
end

function rate_names = rateToReaction(model)
    rate_names = '[';
    reactions = model.get('Reactions');
    nr = length(model.get('Reactions')); 
    for rs = 1:(nr - 1)
        rate_names = strcat(rate_names, '"', getRateName(reactions(rs)), '"', ',' );
    end
    rate_names = strcat(rate_names, '"', getRateName(reactions(nr)), '"]');
end


% isSReac: flag to indicate if this is a surface reactions,
% simply, zero order reactions are not allowed as surface reactions
function [value, factor] = getRateValue(reaction, model, isSReac)
    % check if the reaction got a kinetic law.
    % if not it is a 0-order reaction and
    % take the reaction rate
    % NOTE: STEPS uses mol/l for reaction constants
    if(isempty(reaction.get('KineticLaw')))
        if(isempty(reaction.get('ReactionRate')))
            error('Reaction rate missing: ', reaction)
        end

        name = reaction.get('ReactionRate');

        if(isnumeric(name))
            value = name;
        else
            % get rate from the base workspace
            % this requires that the variable
            % rate_var = evalin('base', name);
            % value = rate_var.get('Value');
            paras = model.get('Parameters');
            para_names = paras.get('Name');
            index = find(strcmp(para_names, name));
            value = paras(index).get('Value');

            % get the value unit
            value_unit = paras(index).get('ValueUnits');
            if(isempty(value_unit))
                error('STEPS for Matlab expects reaction rates with value units.')
            end

            % scale by the unit we got
            factor = 1.0;
            if(strcmp(value_unit,'1/(micromolarity*second)'))
                factor = 1e6;
            elseif(strcmp(value_unit,'1/(nanomolarity*second)'))
                factor = 1e9;
            elseif(strcmp(value_unit, 'micromole/second') & ~isSReac)
                factor = 1e-6; 
            elseif(strcmp(value_unit, 'nanomole/second') & ~isSReac)
                factor = 1e-9;
             elseif(strcmp(value_unit, '1/second'))
                % do nothing - the factor is 1
            else
                error('Unexpected reaction rate unit: ', value_unit);
            end
            value = value * factor;
        end
        return
    end

    % branch on the type of the kinetic law parameter
    % it could either be a double value of a Simbiology Parameter object
    type_str = class(reaction.get('KineticLaw').get('Parameters'));

    if(strcmp(type_str, 'double'))
        value = reaction.get('KineticLaw').get('Parameters');
        if(isempty(value))
            % use the ParameterVariableNames instead
            var_names = reaction.get('KineticLaw').get('ParameterVariableNames');
            if(length(var_names) ~= 1)
                error('Expected only a single parameter for the kinetic law: ', var_names)
            end
            % rate_var = evalin('base', var_names{1});
            % value = rate_var.get('Value') % FIXME: replace with access via model
            paras = model.get('Parameters');
            para_names = paras.get('Name');
            index = find(strcmp(para_names, var_names{1}));
            value = paras(index).get('Value');
 
            % get the value unit
            value_unit = paras(index).get('ValueUnits');
            if(isempty(value_unit))
                error('STEPS for Matlab expects reaction rates with value units.')
            end

            % scale by the unit we got
            factor = 1.0;
            if(strcmp(value_unit,'1/(micromolarity*second)'))
                factor = 1e6;
            elseif(strcmp(value_unit,'1/(nanomolarity*second)'))
                factor = 1e9;
            elseif(strcmp(value_unit, 'micromole/second') & ~isSReac)
                factor = 1e-6; 
            elseif(strcmp(value_unit, 'nanomole/second') & ~isSReac)
                factor = 1e-9;
            elseif(strcmp(value_unit, '1/second'))
                % do nothing - the factor is 1
            else
                error('Unexpected reaction rate unit: ', value_unit);
            end
            value = value * factor;
        end
    elseif(strcmp(type_str, 'SimBiology.Parameter'))
        value = reaction.get('KineticLaw').get('Parameters').get('Value');
    
            % get the value unit
            value_unit = reaction.get('KineticLaw').get('Parameters').get('ValueUnits');
            if(isempty(value_unit))
                error('STEPS for Matlab expects reaction rates with value units.')
            end

            % scale by the unit we got
            factor = 1.0;
            if(strcmp(value_unit,'1/(micromolarity*second)'))
                factor = 1e6;
            elseif(strcmp(value_unit,'1/(nanomolarity*second)'))
                factor = 1e9;
            elseif(strcmp(value_unit, 'micromole/second') & ~isSReac)
                factor = 1e-6; 
            elseif(strcmp(value_unit, 'nanomole/second') & ~isSReac)
                factor = 1e-9;
            elseif(strcmp(value_unit, '1/second'))
                % do nothing - the factor is 1
            else
                error('Unexpected reaction rate unit: ', value_unit);
            end
            value = value * factor;
    else
        error('Parameter is of unknown type: %s.', type_str);
    end
end

function capacity = getCapacityLiter(c)
    % get the compartments capacity
    capacity = c.get('Capacity');
    unit = c.get('CapacityUnits');
        
    % check if the capacity unit is supported
    % NOTE: for reaction rate constants STEPS scales to liter
    factor = 1.0;
    if(strcmp(unit, 'l') || strcmp(unit, 'liter'))
        % do nothing
    elseif(strcmp(unit, 'ml') || strcmp(unit, 'milliliter'))
        factor = 1e3;
    elseif(strcmp(unit, 'mul') || strcmp(unit, 'microliter'))
        factor = 1e6;
    elseif(strcmp(unit, 'nl') || strcmp(unit, 'nanoliter'))
        factor = 1e9;
    elseif(strcmp(unit, ''))
        warning('No capacity unit given, using default of 1 mul.')
        factor = 1e6;
    else
        error('Unsupported capacity unit for compartment.')
    end    
    capacity = capacity * factor;
end

% handle a reaction event
function [reaction_code, activation_code, current_index] = handleReaction(r, rindex, model)
    reaction_code = '';
    activation_code = '';

    eqn = r.get('Reaction');
    stoch = r.get('Stoichiometry');

    % parse the string - only consider non reversible reactions
    if(strfind(eqn, '<->'))
        error('Reversible reactions are not yet supported.')
    end

    % get a cell arreay of the compartment(s)
    comp = getReactionCompartments(r, model);

    % this is a reaction within a compartment,
    % so there should be only one compartment
    % in the list
    if(length(comp) ~= 1)
        error('Reaction is within one compartment. Only one compartment expected: ', comp);
    end

    tmp_code = strcat('r', num2str(rindex), ' = smod.Reac("r', num2str(rindex), '", vsys, lhs = [');
    i = 1;
    first = true;
    % iterate the reactants
    for rs = r.get('Reactants')'
        sname = rs.get('Name');
       
        for j = 1:abs(stoch(i))
            if(first)
                tmp_code = strcat(tmp_code, sname);
                first = false;
            else
                tmp_code = strcat(tmp_code, ', ', sname);
            end
        end

        i = i + 1;
    end

    tmp_code = strcat(tmp_code, '], rhs = [');

    first = true;
    for ps = r.get('Products')'
        sname = ps.get('Name');

        for j = 1:stoch(i)
            if(first)
                tmp_code = strcat(tmp_code , sname);
                first = false;
            else
                tmp_code = strcat(tmp_code, ', ', sname);
            end
        end

        i = i + 1;
    end
    [rate_value, scaling_factor] = getRateValue(r, model, false);
    tmp_code = strcat(tmp_code, '], kcst = ', num2str(rate_value), ' )\n' );
    tmp_code = strcat(tmp_code, 'kcst_si_factor["r', num2str(rindex), '"] = ', num2str(scaling_factor), '\n');     
    compartments = model.get('Compartments');
    for c = compartments'
        activation_code = strcat('sim.setCompReacActive("', c.get('Name'), '", "r', num2str(rindex), '", False)\n');
    end
    activation_code = strcat(activation_code, 'sim.setCompReacActive("',char(comp), '", "r', num2str(rindex), '" , True)\n');
    reaction_code = strcat(reaction_code, tmp_code);
    current_index = rindex;
 
    reaction_code = strcat(reaction_code, strcat('reac_to_comp["r', num2str(rindex), '"] = "', char(comp), '"\n'));

    activation_code = char(activation_code);
end

function [code, activation_code] = getReactions(model)
    code = '';
    activation_code = '';
    
    % only use one index
    index = 1;
    patches = {};
    
    for r = model.get('Reactions')'
        eqn = r.get('Reaction');
        
        % parse the string - only consider non reversible reactions
        if(strfind(eqn, '<->'))
           error('Reversible reactions are not yet supported');
        end

        %  check if we got a reaction or surface reaction
        comp = getReactionCompartments(r, model);
        
        % dispatch to handler functions
        if(length(comp) <= 1)
            [tmp_code, tmp2_code, current_index] = handleReaction(r, index, model);
            code = strcat(code, tmp_code);
            activation_code = [activation_code tmp2_code];
            index = current_index + 1;
        elseif(length(comp) == 2)
            [tmp_code, tmp_patches] = handlePatch(comp, patches);
            code = strcat(code, tmp_code);
            patches = [patches, tmp_patches];
            [tmp_code, tmp2_code] = handleSurfaceReaction(r, index, comp{1}, comp{2}, model, patches);
            code = strcat(code, tmp_code);
            activation_code = [activation_code tmp2_code];
            index = index + 1;
        else
            error('Unexpected number of compartments.')
        end
    end
end

% extract the events
function event_code = getEvents(model)
    % get all events from the model
    events = model.get('Events');

    event_str = {};
    counter = 0;
    for i = 1:length(events)
        % for each event get the trigger and the function
   
        % check that they are only time based,
        % since those are the only ones we support
       
        % get the event functions
        events(i).get('EventFcns');
        fcns = char(events(i).get('EventFcns'));

        % get the event triggers
        trgs = char(events(i).get('Trigger'));
        trgs = trgs(~isspace(trgs));

        % all further checks on the python side
        event_str= [event_str, struct(strcat('trigger_', num2str(counter)), trgs, strcat('event_', num2str(counter)), fcns)];

        counter = counter + 1;
    end

    % JSON encode
    event_code = jsonencode(event_str);

    % rest of the event handling such as building
    % an event queue and ensuring that we handle them at the
    % requested points in time is on the Python side
end

function [model_code, set_initial_conditions] = extractModel(model)
    % preamble
    code = 'import steps.model as smod\n';
    code = strcat(code, 'import steps.geom as geom\n');
    code = strcat(code, 'model = smod.Model()\n');
    code = strcat(code, 'vsys = smod.Volsys("vsys", model)\n');
    code = strcat(code, 'geometry = swm.Geom()\n');
    code = strcat(code, 'ssys = smod.Surfsys("ssys", model)\n');

    % simple lists to get easy access to all species and compartments
    code = strcat(code, 'spec_names = []\n');
    code = strcat(code, 'comp_names = []\n');
    code = strcat(code, 'reac_to_comp = {}\n');
    code = strcat(code, 'sreac_to_patch = {}\n');
    code = strcat(code, 'kcst_si_factor = {}\n');
    code = strcat(code, 'comp_si_factor = {}\n');

    % get compartments - not the nested ones
    compartments = model.get('Compartments');
    
    % get a flat list of all the nested compartments
    set_compartments = '';
    set_species = '';
    set_initial_conditions = '';
    set_reactions = '';
    list_of_compartments = zeros(2, 1);
    first_call = true;
    compartment_capacities = [];
    for compartment = compartments'
        [c, a, b, v] = getCompartments(compartment, model, first_call);
        set_compartments = strcat(set_compartments, c);
        set_species = strcat(set_species, a);
        set_initial_conditions = strcat(set_initial_conditions, b);
        first_call = false;
    end
    
    % create the reactions
    [tmp_code, tmp2_code] = getReactions(model);
    set_reactions = strvcat(set_reactions, tmp_code);

    model_code = strcat(code, set_compartments, set_species, set_reactions);

    % attach the activation of reactions to the initial conditions
    set_initial_conditions = strcat(set_initial_conditions, tmp2_code);
end
% done extract model
