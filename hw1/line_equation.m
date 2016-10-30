function eqn = line_equation(pt1, pt2, method)
    if nargin < 3
        method = 1;
    end
    
	if method == 1
% 	    pt1 = pt1 / sum(pt1);
% 	    pt2 = pt2 / sum(pt1);
	    if length(pt1) < 3
	        pt1 = [pt1; 1];
	        pt2 = [pt2; 1];
	    end
	    
	    eqn = cross(pt1, pt2);
	else
	    x_diff = pt2(1) - pt1(1);
	    if x_diff == 0
	    	eqn = [1; 0; -1 * pt1(1)];
	        return;
	    end

	    % slope
	    m = (pt2(2) - pt1(2)) / x_diff;

	    a = m;
	    b = -1;
	    c = pt1(2) - m * pt1(1);
	    eqn = [a; b; c];
	    % eqn = eqn / norm(eqn, 1);
	end
end
