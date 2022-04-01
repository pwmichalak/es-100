% obtain data from excel sheet
function channel = get_data(file, sheet)

    % obtain excel DST data in a MATLAB-friendly format
    % file: data file name
    % sheet: name of the sheet in the excel file
    % 
    % Conventions used:
    % y: voltage (output) data
    % u: current (input) data
    % t: time range
    % channel: the dataset
    
    [C,columns] = xlsread(file,sheet); % file may need to be .xlsx

    % clean up the column names
    columns = strrep(columns, '(', '_');
    columns = strrep(columns, ')', '');
    columns = strrep(columns, ' ', '');
    columns = strrep(columns, '/', '_');
    
    channel = array2table(C, 'VariableNames', columns);
end