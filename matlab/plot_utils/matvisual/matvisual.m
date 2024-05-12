% +------------------------------------------------------+
% |   Matrix Visualization with MATLAB Implementation    | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        12/02/18 | 
% +------------------------------------------------------+
% 
% function: matvisual(A, varargin)
%
% Input:
% A - m-by-n-by-p matrix to be visualized
% varargin - type 'annotation' in the place of varargin if one wants to
%            annotate the matrix plot (x-label, y-label, etc.)

function matvisual(A, varargin)

% input validation
validateattributes(A, {'single', 'double', 'logical'}, ...
                      {'3d', 'nonempty', 'real'}, ...
                      '', 'A', 1)

% determine the matrix size
[M, N, P] = size(A);

% loop through the matrix pages
for p = 1:P
    % get prepared for page-by-page visualization
    if P > 1, subplot(1, P, p), end 
    
    % visualize the matrix
    himg = imagesc(A(:, :, p));
    grid on
      
    % values labeling
    for m = 1:M
        for n = 1:N
            text(n, m, num2str(A(m, n, p), 3), ...
                'FontName', 'Times New Roman', ...
                'FontSize', round(6 + 50./sqrt(M.*N)), ...
                'HorizontalAlignment', 'center', ...
                'Rotation', 45)
        end
    end
    
    % annotation
    if strcmp(varargin, 'annotation')
        % x-label, y-label, x-ticks, y-ticks, title
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
        xlabel('Column number')
        ylabel('Row number')
        if P > 1, title(['Matrix page ' num2str(p)]), end
        if M <= 50, set(gca, 'YTick', 1:M), end
        if N <= 50, set(gca, 'XTick', 1:N), end
    end
    
    % set the datatip UpdateFcn
    cursorMode = datacursormode(gcf);
    set(cursorMode, 'UpdateFcn', {@datatiptxt, himg})
end

end

function text_to_display = datatiptxt(~, hDatatip, himg)

% determine the current datatip position
pos = get(hDatatip, 'Position');

% form the datatip label
text_to_display = {['Row: ' num2str(pos(2))], ...
                   ['Column: ' num2str(pos(1))], ...
                   ['Value: ' num2str(himg.CData(pos(2), pos(1)))]};

end