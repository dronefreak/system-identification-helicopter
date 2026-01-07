function reportFile = generate_report(results, varargin)
% GENERATE_REPORT - Generate comprehensive optimization report
%
% This function creates automated reports from optimization results including
% plots, statistics, and formatted summaries in multiple formats.
%
% Usage:
%    generate_report(results)
%    generate_report(results, 'format', 'pdf')
%    reportFile = generate_report(results, 'outputDir', 'reports/')
%
% Inputs:
%    results      - Results structure from run_experiment or IWO output
%
% Optional Parameters (Name-Value pairs):
%    'format'       - 'html', 'pdf', 'markdown', or 'all' (default: 'html')
%    'outputDir'    - Output directory (default: 'reports/')
%    'filename'     - Base filename (default: auto from results)
%    'include'      - Cell array: {'plots', 'stats', 'params', 'metadata'}
%                     (default: all)
%    'plotFormat'   - 'png', 'pdf', 'fig' (default: 'png')
%    'openReport'   - Open report after generation (default: true)
%
% Outputs:
%    reportFile     - Path to generated report
%
% Report Contents:
%    - Executive summary
%    - Convergence plots
%    - Statistical analysis
%    - Parameter values
%    - Metadata and configuration
%    - System information
%
% Examples:
%    % Generate HTML report
%    generate_report(results);
%
%    % Generate PDF report
%    generate_report(results, 'format', 'pdf');
%
%    % Custom output with specific sections
%    generate_report(results, ...
%                    'outputDir', 'reports/batch1/', ...
%                    'include', {'plots', 'stats'});
%
% See also: run_experiment, plot_convergence, plot_comparison
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Inputs
    p = inputParser;
    addRequired(p, 'results', @isstruct);
    addParameter(p, 'format', 'html', @(x) any(validatestring(x, {'html', 'pdf', 'markdown', 'all'})));
    addParameter(p, 'outputDir', 'reports/', @ischar);
    addParameter(p, 'filename', '', @ischar);
    addParameter(p, 'include', {'plots', 'stats', 'params', 'metadata'}, @iscell);
    addParameter(p, 'plotFormat', 'png', @(x) any(validatestring(x, {'png', 'pdf', 'fig'})));
    addParameter(p, 'openReport', true, @islogical);
    parse(p, results, varargin{:});

    results = p.Results.results;
    format = p.Results.format;
    outputDir = p.Results.outputDir;
    filename = p.Results.filename;
    includeSections = p.Results.include;
    plotFormat = p.Results.plotFormat;
    openReport = p.Results.openReport;

    fprintf('=======================================================\n');
    fprintf('  Generating Optimization Report\n');
    fprintf('=======================================================\n\n');

    %% Prepare Output Directory
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(scriptPath));

    if ~isempty(outputDir) && outputDir(1) ~= '/' && ~(length(outputDir) > 1 && outputDir(2) == ':')
        outputDir = fullfile(projectRoot, outputDir);
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
        fprintf('Created output directory: %s\n', outputDir);
    end

    %% Determine Filename
    if isempty(filename)
        if isfield(results, 'metadata') && isfield(results.metadata, 'experimentName')
            filename = results.metadata.experimentName;
        else
            filename = sprintf('report_%s', datestr(now, 'yyyymmdd_HHMMSS'));
        end
    end

    %% Extract Data
    fprintf('Extracting data...\n');

    reportData = struct();

    % Results data
    if isfield(results, 'bestSolution')
        reportData.bestSolution = results.bestSolution;
    elseif isfield(results, 'BestSol')
        reportData.bestSolution = results.BestSol;
    end

    if isfield(results, 'bestCosts')
        reportData.bestCosts = results.bestCosts;
    elseif isfield(results, 'BestCosts')
        reportData.bestCosts = results.BestCosts;
    end

    if isfield(results, 'metadata')
        reportData.metadata = results.metadata;
    end

    if isfield(results, 'config')
        reportData.config = results.config;
    end

    if isfield(results, 'timing')
        reportData.timing = results.timing;
    end

    %% Generate Plots
    figureDir = fullfile(outputDir, 'figures');
    if any(strcmp(includeSections, 'plots'))
        if ~exist(figureDir, 'dir')
            mkdir(figureDir);
        end

        fprintf('Generating plots...\n');

        % Convergence plot
        if isfield(reportData, 'bestCosts')
            convergenceFig = plot_convergence(reportData.bestCosts, ...
                                             'title', 'Optimization Convergence');
            convergencePath = fullfile(figureDir, ['convergence.' plotFormat]);
            saveas(convergenceFig, convergencePath);
            close(convergenceFig);
            reportData.convergencePlot = convergencePath;
            fprintf('  ✓ Convergence plot saved\n');
        end
    end

    %% Generate Report Based on Format
    fprintf('\nGenerating %s report...\n', upper(format));

    if strcmp(format, 'all')
        formats = {'html', 'markdown'};
    else
        formats = {format};
    end

    reportFiles = {};
    for i = 1:length(formats)
        fmt = formats{i};

        switch fmt
            case 'html'
                reportFile = generate_html_report(reportData, outputDir, filename, includeSections);
                reportFiles{end+1} = reportFile; %#ok<AGROW>

            case 'markdown'
                reportFile = generate_markdown_report(reportData, outputDir, filename, includeSections);
                reportFiles{end+1} = reportFile; %#ok<AGROW>

            case 'pdf'
                fprintf('  PDF generation requires publish() or external tool\n');
                fprintf('  Generating Markdown instead, convert with pandoc:\n');
                fprintf('  pandoc report.md -o report.pdf\n');
                reportFile = generate_markdown_report(reportData, outputDir, filename, includeSections);
                reportFiles{end+1} = reportFile; %#ok<AGROW>
        end
    end

    fprintf('\n=======================================================\n');
    fprintf('Report generated successfully!\n');
    for i = 1:length(reportFiles)
        fprintf('  %s\n', reportFiles{i});
    end
    fprintf('=======================================================\n\n');

    %% Open Report
    if openReport && ~isempty(reportFiles)
        try
            if ispc
                system(['start "" "' reportFiles{1} '"']);
            elseif ismac
                system(['open "' reportFiles{1} '"']);
            elseif isunix
                system(['xdg-open "' reportFiles{1} '" &']);
            end
        catch
            fprintf('Could not auto-open report. Please open manually.\n');
        end
    end

    if ~isempty(reportFiles)
        reportFile = reportFiles{1};
    else
        reportFile = '';
    end
end

function reportFile = generate_html_report(reportData, outputDir, filename, includeSections)
    % Generate HTML report
    reportFile = fullfile(outputDir, [filename '.html']);

    fid = fopen(reportFile, 'w');
    if fid == -1
        error('GenerateReport:FileError', 'Could not create report file');
    end

    % HTML header
    fprintf(fid, '<!DOCTYPE html>\n');
    fprintf(fid, '<html>\n<head>\n');
    fprintf(fid, '<title>Optimization Report: %s</title>\n', filename);
    fprintf(fid, '<style>\n');
    fprintf(fid, '  body { font-family: Arial, sans-serif; margin: 40px; }\n');
    fprintf(fid, '  h1 { color: #333; }\n');
    fprintf(fid, '  h2 { color: #666; border-bottom: 2px solid #ddd; padding-bottom: 5px; }\n');
    fprintf(fid, '  table { border-collapse: collapse; width: 100%%; margin: 20px 0; }\n');
    fprintf(fid, '  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n');
    fprintf(fid, '  th { background-color: #f2f2f2; }\n');
    fprintf(fid, '  .metric { font-weight: bold; color: #0066cc; }\n');
    fprintf(fid, '  img { max-width: 100%%; height: auto; margin: 20px 0; }\n');
    fprintf(fid, '</style>\n');
    fprintf(fid, '</head>\n<body>\n');

    % Title
    fprintf(fid, '<h1>Optimization Report</h1>\n');
    if isfield(reportData, 'metadata') && isfield(reportData.metadata, 'experimentName')
        fprintf(fid, '<p><strong>Experiment:</strong> %s</p>\n', reportData.metadata.experimentName);
    end
    fprintf(fid, '<p><strong>Generated:</strong> %s</p>\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    % Executive Summary
    fprintf(fid, '<h2>Executive Summary</h2>\n');
    if isfield(reportData, 'metadata')
        fprintf(fid, '<table>\n');
        if isfield(reportData.metadata, 'finalCost')
            fprintf(fid, '<tr><td>Final Cost</td><td class="metric">%.6f</td></tr>\n', reportData.metadata.finalCost);
        end
        if isfield(reportData.metadata, 'initialCost')
            fprintf(fid, '<tr><td>Initial Cost</td><td>%.6f</td></tr>\n', reportData.metadata.initialCost);
        end
        if isfield(reportData.metadata, 'improvement')
            fprintf(fid, '<tr><td>Improvement</td><td class="metric">%.2f%%</td></tr>\n', reportData.metadata.improvement);
        end
        if isfield(reportData.metadata, 'executionTime')
            fprintf(fid, '<tr><td>Execution Time</td><td>%.2f minutes</td></tr>\n', reportData.metadata.executionTime / 60);
        end
        fprintf(fid, '</table>\n');
    end

    % Plots
    if any(strcmp(includeSections, 'plots')) && isfield(reportData, 'convergencePlot')
        fprintf(fid, '<h2>Convergence Plot</h2>\n');
        [~, relPath] = fileparts(reportData.convergencePlot);
        fprintf(fid, '<img src="figures/%s" alt="Convergence Plot">\n', [relPath '.png']);
    end

    % Statistics
    if any(strcmp(includeSections, 'stats')) && isfield(reportData, 'bestCosts')
        fprintf(fid, '<h2>Statistics</h2>\n');
        fprintf(fid, '<table>\n');
        fprintf(fid, '<tr><td>Total Iterations</td><td>%d</td></tr>\n', length(reportData.bestCosts));
        fprintf(fid, '<tr><td>Best Cost</td><td>%.6f</td></tr>\n', min(reportData.bestCosts));
        fprintf(fid, '<tr><td>Worst Cost</td><td>%.6f</td></tr>\n', max(reportData.bestCosts));
        fprintf(fid, '</table>\n');
    end

    % Parameters
    if any(strcmp(includeSections, 'params')) && isfield(reportData, 'bestSolution')
        fprintf(fid, '<h2>Best Parameters</h2>\n');
        fprintf(fid, '<table>\n');
        fprintf(fid, '<tr><th>Index</th><th>Value</th></tr>\n');
        params = reportData.bestSolution.Position;
        for i = 1:min(40, length(params))
            fprintf(fid, '<tr><td>%d</td><td>%.6f</td></tr>\n', i, params(i));
        end
        fprintf(fid, '</table>\n');
    end

    % Metadata
    if any(strcmp(includeSections, 'metadata')) && isfield(reportData, 'metadata')
        fprintf(fid, '<h2>Metadata</h2>\n');
        fprintf(fid, '<table>\n');
        meta = reportData.metadata;
        if isfield(meta, 'platform')
            fprintf(fid, '<tr><td>Platform</td><td>%s</td></tr>\n', meta.platform);
        end
        if isfield(meta, 'matlabRelease')
            fprintf(fid, '<tr><td>MATLAB</td><td>%s</td></tr>\n', meta.matlabRelease);
        end
        fprintf(fid, '</table>\n');
    end

    % Footer
    fprintf(fid, '<hr>\n');
    fprintf(fid, '<p><small>Generated by System Identification Toolkit</small></p>\n');
    fprintf(fid, '</body>\n</html>\n');

    fclose(fid);
end

function reportFile = generate_markdown_report(reportData, outputDir, filename, includeSections)
    % Generate Markdown report
    reportFile = fullfile(outputDir, [filename '.md']);

    fid = fopen(reportFile, 'w');
    if fid == -1
        error('GenerateReport:FileError', 'Could not create report file');
    end

    % Title
    fprintf(fid, '# Optimization Report\n\n');
    if isfield(reportData, 'metadata') && isfield(reportData.metadata, 'experimentName')
        fprintf(fid, '**Experiment:** %s\n\n', reportData.metadata.experimentName);
    end
    fprintf(fid, '**Generated:** %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    % Executive Summary
    fprintf(fid, '## Executive Summary\n\n');
    if isfield(reportData, 'metadata')
        if isfield(reportData.metadata, 'finalCost')
            fprintf(fid, '- **Final Cost:** %.6f\n', reportData.metadata.finalCost);
        end
        if isfield(reportData.metadata, 'initialCost')
            fprintf(fid, '- **Initial Cost:** %.6f\n', reportData.metadata.initialCost);
        end
        if isfield(reportData.metadata, 'improvement')
            fprintf(fid, '- **Improvement:** %.2f%%\n', reportData.metadata.improvement);
        end
        if isfield(reportData.metadata, 'executionTime')
            fprintf(fid, '- **Execution Time:** %.2f minutes\n', reportData.metadata.executionTime / 60);
        end
        fprintf(fid, '\n');
    end

    % Plots
    if any(strcmp(includeSections, 'plots')) && isfield(reportData, 'convergencePlot')
        fprintf(fid, '## Convergence Plot\n\n');
        [~, relPath] = fileparts(reportData.convergencePlot);
        fprintf(fid, '![Convergence](figures/%s)\n\n', [relPath '.png']);
    end

    % Statistics
    if any(strcmp(includeSections, 'stats')) && isfield(reportData, 'bestCosts')
        fprintf(fid, '## Statistics\n\n');
        fprintf(fid, '| Metric | Value |\n');
        fprintf(fid, '|--------|-------|\n');
        fprintf(fid, '| Total Iterations | %d |\n', length(reportData.bestCosts));
        fprintf(fid, '| Best Cost | %.6f |\n', min(reportData.bestCosts));
        fprintf(fid, '| Worst Cost | %.6f |\n', max(reportData.bestCosts));
        fprintf(fid, '\n');
    end

    % Parameters
    if any(strcmp(includeSections, 'params')) && isfield(reportData, 'bestSolution')
        fprintf(fid, '## Best Parameters\n\n');
        fprintf(fid, '| Index | Value |\n');
        fprintf(fid, '|-------|-------|\n');
        params = reportData.bestSolution.Position;
        for i = 1:min(40, length(params))
            fprintf(fid, '| %d | %.6f |\n', i, params(i));
        end
        if length(params) > 40
            fprintf(fid, '\n*... and %d more parameters*\n', length(params) - 40);
        end
        fprintf(fid, '\n');
    end

    % Metadata
    if any(strcmp(includeSections, 'metadata')) && isfield(reportData, 'metadata')
        fprintf(fid, '## Metadata\n\n');
        fprintf(fid, '| Property | Value |\n');
        fprintf(fid, '|----------|-------|\n');
        meta = reportData.metadata;
        if isfield(meta, 'platform')
            fprintf(fid, '| Platform | %s |\n', meta.platform);
        end
        if isfield(meta, 'matlabRelease')
            fprintf(fid, '| MATLAB | %s |\n', meta.matlabRelease);
        end
        if isfield(reportData, 'timing')
            fprintf(fid, '| Execution Time | %.2f seconds |\n', reportData.timing);
        end
        fprintf(fid, '\n');
    end

    % Configuration
    if any(strcmp(includeSections, 'metadata')) && isfield(reportData, 'config')
        fprintf(fid, '## Configuration\n\n');
        fprintf(fid, '```\n');
        config = reportData.config;
        fields = fieldnames(config);
        for i = 1:length(fields)
            value = config.(fields{i});
            if isnumeric(value)
                if isscalar(value)
                    fprintf(fid, '%s: %g\n', fields{i}, value);
                else
                    fprintf(fid, '%s: [%s]\n', fields{i}, num2str(value, '%.2f '));
                end
            elseif ischar(value)
                fprintf(fid, '%s: %s\n', fields{i}, value);
            elseif islogical(value)
                fprintf(fid, '%s: %s\n', fields{i}, mat2str(value));
            end
        end
        fprintf(fid, '```\n\n');
    end

    % Footer
    fprintf(fid, '---\n\n');
    fprintf(fid, '*Generated by System Identification Toolkit*\n');

    fclose(fid);
end
