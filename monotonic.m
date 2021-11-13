% The following function gets the time interval, y0, yref (y_target), the coefficient alpha,
% the control c, the matrices A and B as well as the mesh size N, and updates the control c
% using the monotonic algorithm and returns back the updated control!

function [control] = monotonic(interval_start, interval_end, y0, yref, alpha, c, A, B, N)

  t = linspace(interval_start, interval_end, N);
  dt = t(2)-t(1);

% Initializations of y and p
  y = zeros(2,N);
  y(:,1) = y0;
  for n=1:N-1
    y(:,n+1) = expm(i*A*dt)*expm(i*B*c(n)*dt)*y(:,n);
  end
  
  p = zeros(2,N);
  p(:,N) = -yref;
  for n=N:-1:2
    p(:,n-1) = expm(-i*B*c(n-1)*dt)*expm(-i*A*dt)*p(:,n);
  end
% End of the initialization

  lambda = 1;

  for n=1:N-1
    c(n) = solvemono(c, n, alpha, dt, B, A, y, p, lambda);
    y(:,n+1)=expm(i*A*dt)*expm(i*B*c(n)*dt)*y(:,n);
  end
  
  control = c;
  
end



% The following function solves the equation:
% x^2-c(n)^2 + 2/(alpha*dt) * (real(((expm(i*(x-c(n))*B*dt)-eye(size(A)))*y(:,n))'*p(:,n))) + lambda*(x-c(n))^2 = 0
% We use the Newton method [x_new = x_old - f(x_old)/df(x_old)] for the following f:
% f(x) = x^2-c(n)^2 + 2/(alpha*dt) * (real(((expm(i*(x-c(n))*B*dt)-eye(size(A)))*y(:,n))'*p(:,n))) + lambda*(x-c(n))^2 
% But this may fail, because x = c(n) is a trivial root of the equation!
% To avoid such problem, we let the factor x-c(n) out!
% f(x) = x+c(n) + 2/(alpha*dt*(x-c(n))) * (real(((expm(i*(x-c(n))*B*dt)-eye(size(A)))*y(:,n))'*p(:,n))) + lambda*(x-c(n))

function x = solvemono(c, n, alpha, dt, B, A, y, p, lambda)
  
  err = 1;
  
  if n == 1;
    x = c(n) + 1e-12;
  else
    x = c(n-1); %Two natural ways to initialize: x=c^k(n) or x = c^{k+1}(n-1).
  end
  
  while err > (10e-5)
    x_old = x;
    x = x - (x+c(n) + 2/(alpha*dt*(x-c(n))) * (real(((expm(i*(x-c(n))*B*dt)-eye(size(A)))*y(:,n))'*p(:,n))) + lambda*(x-c(n))) ...
            /(1 - alpha*dt*2/(alpha*dt*(x-c(n)))^2 * (real(((expm(i*(x-c(n))*B*dt)-eye(size(A)))*y(:,n))'*p(:,n))) ...
	            + 2/(alpha*dt*(x-c(n))) * (real((((i*B*dt)*expm(i*(x-c(n))*B*dt))*y(:,n))'*p(:,n))) + lambda);
    err = abs(x-x_old);
  end

end
