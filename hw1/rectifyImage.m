function [rectI, H, H_aff] = rectifyImage(filename, debug)

    if nargin == 1
        debug = 0;
    end
    
    [d, fname] = fileparts(filename);
    mat_fpath = fullfile(strrep(d, 'images', 'data'), [fname, '.mat']);
    if debug == 0
        load(mat_fpath);
    end

    img = imread(filename);
    
    if debug == 1
        % annotate parallel points
        [x, y] = get_points(img, 4);
    end
    
    [rect_affine, H_aff] = rectify_affine(img, x, y, filename);
    
    imshow(rect_affine);
    waitforbuttonpress; close;
    aff_fname = [filename(1:end-4), '-affine.jpg'];
    imwrite(rect_affine, aff_fname)
    
    if debug == 1
        % annotate perpendicular points
        [px, py] = get_points(rect_affine, 12);
        save(mat_fpath, 'x', 'y', 'px', 'py');
    end
    
    [rectI, H] = rectify_euclidean(rect_affine, px, py, filename);
    H = inv(H) * H_aff;
    
    imshow(rectI);
    waitforbuttonpress; close;
    aff_fname = [filename(1:end-4), '-rectified.jpg'];
    imwrite(rectI, aff_fname)
end

%% affine rectification on img
function [rectI, H] = rectify_affine(img, x, y, fname)

    line_a = line_equation([x(1); y(1)], [x(2); y(2)], 1);
    line_b = line_equation([x(3); y(3)], [x(4); y(4)], 1);
    line_c = line_equation([x(1); y(1)], [x(4); y(4)], 1);
    line_d = line_equation([x(2); y(2)], [x(3); y(3)], 1);
    
    
    inf_p1 = line_intersect(line_a, line_b);
    inf_p2 = line_intersect(line_c, line_d);
    
    line_inf = line_equation(inf_p1, inf_p2, 0);

    w = size(img, 2);
    llel_fname = [fname(1:end-4), '-llel.jpg'];
%     range = [min([-w, inf_p1(1), inf_p2(1)]), max([w, inf_p1(1), inf_p2(1)])];
    range = [-w, w*2];
    visualize_lines(img, {line_a, line_b, line_c, line_d, line_inf}, llel_fname, range)
        
    H = eye(3);
    H(3,:) = line_inf;
    rectI = applyH(img, H);
end

%% Euclidean rectification
function [rectI, H] = rectify_euclidean(img, x, y, fname)
    line_a = line_equation([x(1); y(1)], [x(2); y(2)]);
    line_b = line_equation([x(3); y(3)], [x(4); y(4)]);
    line_c = line_equation([x(5); y(5)], [x(6); y(6)]);
    line_d = line_equation([x(7); y(7)], [x(8); y(8)]);
    line_e = line_equation([x(9); y(9)], [x(10); y(10)]);
    line_f = line_equation([x(11); y(11)], [x(12); y(12)]);
    
    
    w = size(img, 2);
    perp_fname = [fname(1:end-4), '-perp.jpg'];
    visualize_lines(img, {line_a, line_b, line_c, line_d}, perp_fname, [0, w])
    
    A = [line_a(1) * line_b(1), line_a(1) * line_b(2) + line_a(2) * line_b(1);
         line_c(1) * line_d(1), line_c(1) * line_d(2) + line_c(2) * line_d(1);
         ];
     
    b = [-line_a(2) * line_b(2);
         -line_c(2) * line_d(2);
        ];
    
    X = A \ b;
    
    L = [X(1), X(2);
         X(2), 1];
     
    [U, S, ~] = svd(L);
    
    U = U * sqrt(S) * U';
    
    H = eye(3);
    H(1:2, 1:2) = U;
    
    rectI = applyH(img, inv(H));
    
    fprintf('\nBefore: %.7f, After: %8.7E', ...
        cos_transform(eye(3), line_a, line_b), cos_transform(H, line_a, line_b));
    fprintf('\nBefore: %.7f, After: %8.7E', ...
        cos_transform(eye(3), line_c, line_d), cos_transform(H, line_c, line_d));
    fprintf('\nBefore: %.7f, After: %8.7E\n\n', ...
        cos_transform(eye(3), line_e, line_f), cos_transform(H, line_e, line_f));
end

%% get points from image
function [x, y] = get_points(img, n)
    imshow(img);
    hold on;
    
    [x, y] = ginput(n);
    hold off;
    close;
end


%% Find cosine angle
function [cos_theta] = cos_transform(H, l1, l2)
    C = eye(3);
    C(3,3) = 0;
    
    C = H * C * H';
    
    num = l1' * C * l2;
    den = sqrt(l1' * C * l1) * sqrt(l2' * C * l2);
    
    cos_theta = num / den;
end

%% visualize all lines on image
function visualize_lines(img, lines, fig_fname, range)
    figure;
    imshow(img);
    hold on;
    
    for i = 1 : length(lines)
        plot_line(lines{i}, range);
    end
    hold off;
    saveas(gcf, fig_fname);
    waitforbuttonpress;
    close;
end


%% line_intersect: Finds the intersection point of two lines
function [point] = line_intersect(l1, l2)
    point = cross(l1, l2);
    point = point / point(3);
end

%% plot_line: Plot a line by its equation
function plot_line(l, range)
    % y = -ax/b -c/b
    f = @(x) (-l(1) / l(2) * x - l(3) / l(2));
    fplot(f, range);
end
