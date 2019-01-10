classdef MyPolyCurve
   properties
       poly_rank
       coeff_matrix
   end
   methods
       function obj = MyPolyCurve(matrix)
           obj.coeff_matrix = matrix;
           obj.poly_rank = size(matrix, 1);
       end
       
      function y = getY(obj, x)
         x_vector = zeros(1, obj.poly_rank + 1);
         for rank = 0 : obj.poly_rank
             x_vector(rank) = x^rank;
         end
         y_poly = x_vector * obj.coeff_matrix;
         y = roots(y_poly);
      end
      
      function x = getX(obj,y)
          y_vector = zeros(obj.poly_rank + 1, 1);
         for rank = 0 : obj.poly_rank
             y_vector(rank) = y^rank;
         end
         x_poly = (obj.coeff_matrix * y_vector).';
         x = roots(x_poly);
      end
      
      function value = getValue(obj, x, y)
          x = x(:).'; % turn into a row vector
          y = y(:); % turn into a column vector
          value = x * obj.coeff_matrix * y;
      end
   end
end