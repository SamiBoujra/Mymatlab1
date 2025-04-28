function rhoDE_global = two_step_wiener_min(M_new, S_noisy, Niter, opts)
% TWO_STEP_WIENER_MIN   Tikhonov‐Wiener reconstruction with PSD projection
%   M_new    : N x NE x Ntau array of forward operators
%   S_noisy  : N x Ntau measured data
%   Niter    : number of iterations per slice
%   opts.gamma_range = [gamma_min, gamma_max]  (default [0,1e-5])
%   opts.nGamma      = number of gamma samples   (default 50)
%   opts.projPSD     = true/false to enforce PSD+unit trace
%   opts.jitter_frac = fraction of ||A||_F for jitter if chol fails

  % defaults
  if ~isfield(opts,'gamma_range'), opts.gamma_range = [0 1e-5]; end
  if ~isfield(opts,'nGamma'),      opts.nGamma      = 50;     end
  if ~isfield(opts,'projPSD'),     opts.projPSD     = true;   end
  if ~isfield(opts,'jitter_frac'), opts.jitter_frac = 1e-8;   end

  [N, Ntau] = size(S_noisy);
  NE        = size(M_new,2);
  rhoDE_global = zeros(NE, Ntau);

  % build gamma list
  gamma_list = logspace(opts.gamma_range(1), opts.gamma_range(2), opts.nGamma);

  mid = ceil(Ntau/2);
  for iDE = 1:mid
    M      = M_new(:,:,iDE);
    S      = S_noisy(:,iDE);
    x_prev = zeros(NE,1);

    bestRes = Inf;
    bestRho = x_prev;

    for iter = 1:Niter
      Rnorm = zeros(opts.nGamma,1);
      Xnorm = zeros(opts.nGamma,1);

      % compute residuals for each gamma
      for j = 1:opts.nGamma
        g = gamma_list(j);
        A = M'*M + g*eye(NE);
        b = M'*S + g*x_prev;

        [R,p] = chol(A);
        if p>0
          A = A + opts.jitter_frac*norm(A,'fro')*eye(NE);
          R = chol(A);
        end

        xg = R \ (R'\ b);
        Rnorm(j) = norm(M*xg - S);
        Xnorm(j) = norm(xg);
      end

      % pick gamma at max curvature of L-curve
      k = lcurve_corner(Rnorm, Xnorm);
      alpha = gamma_list(k);

      % final solve with alpha
      A = M'*M + alpha*eye(NE);
      b = M'*S + alpha*x_prev;
      [R,p] = chol(A);
      if p>0
        A = A + opts.jitter_frac*norm(A,'fro')*eye(NE);
        R = chol(A);
      end
      x_new = R \ (R'\ b);
       %–– Compute residuals for each gamma ––
      Rnorm = zeros(opts.nGamma,1);
      Xnorm = zeros(opts.nGamma,1);
      for j = 1:opts.nGamma
        g = gamma_list(j);
        A = M'*M + g*eye(NE);
        b = M'*S + g*x_prev;
        [R,p] = chol(A);
        if p>0
          A = A + opts.jitter_frac*norm(A,'fro')*eye(NE);
          R = chol(A);
        end
        xg = R \ (R'\ b);
        Rnorm(j) = norm(M*xg - S);
        Xnorm(j) = norm(xg);
      end

      tau =1.01
      delta=1e-5
      %–– DEBUG: plot residual vs gamma to check discrepancy principle ––
      figure(101); clf;
      semilogy(gamma_list, Rnorm, '-o', 'MarkerSize', 6);
      hold on;
      yline(tau*delta, 'r--', 'LineWidth', 1.5);
      xlabel('gamma');
      ylabel('||M x(\gamma) - S||');
      title(sprintf('Slice %d: Residual vs gamma', iDE));
      drawnow;

      % PSD + unit trace projection
      if opts.projPSD
        Xmat = reshape((x_new + x_new')/2, NE, NE);
        [V,D] = eig(Xmat);
        D = diag(max(diag(D),0));
        Xpsd = V * D * V';
        Xpsd = Xpsd / trace(Xpsd);
        x_new = Xpsd(:);
      end

      res = norm(M*x_new - S)/norm(S);
      if res < bestRes
        bestRes = res;
        bestRho = x_new;
      end

      if norm(x_new - x_prev)/max(norm(x_prev),eps) < 1e-6
        break
      end

      x_prev = x_new;
    end

    rhoDE_global(:,iDE) = bestRho;
    fprintf('Slice %d/%d  best residual = %.2e\n', iDE, Ntau, bestRes);
  end

  % mirror the second half
  for iDE = mid+1:Ntau
    rhoDE_global(:,iDE) = rhoDE_global(:, Ntau - iDE + 1);
  end
end


function k = lcurve_corner(rnorm, xnorm)
  % approximate L-curve corner by max discrete curvature
  t = log(rnorm);
  s = log(xnorm);
  dt = diff(t);
  ds = diff(s);
  d2t = diff(dt);
  d2s = diff(ds);
  num = abs(d2t.*ds(1:end-1) - dt(1:end-1).*d2s);
  den = (dt(1:end-1).^2 + ds(1:end-1).^2).^(3/2) + eps;
  [~, idx] = max(num./den);
  k = idx + 1;
end
