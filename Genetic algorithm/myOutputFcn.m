function [state, options, optchanged] = myOutputFcn(options, state, flag)
    % Persistent variable to store the history of best objective values
    persistent history;
    
    if isempty(history)
        history = [];
    end
    
    % Append the best objective value of the current generation to history
    if strcmp(flag, 'iter')
        history = [history; min(state.Score)];
    end
    
    % Update the output structure
    optchanged = false;
    
    % Store the history in the base workspace to access it after optimization
    assignin('base', 'GA_Convergence_History', history);
end
