function [rec_v, rec_sp] = run_LIF( Iinj, stop)
  % Simulate the LIF dynamics with external input current
  % Set parameters
  V_th = 15
  V_reset = 0
  tau_m = 10
  g_L = 1000
  V_init = 0
  E_L = 0
  dt = 0.1
  range_t = 100
  Lt = range_t
  tref = 0.1
  % Initialize voltage
  v = zeros(Lt, 1);
  v(1) = V_init
  % Set current time course
  Iinj = Iinj * ones(Lt, 1)
  % If current pulse, set beginning and end to 0
  if stop == 1
    Iinj(1:floor(length(Iinj) / 2) - 1000) = 0;
    Iinj(floor(length(Iinj) / 2) + 1000:end) = 0;
  end
  % Loop over time
  rec_spikes = []; % record spike times
  tr = 0; % the count for refractory duration
  for it = 1:Lt - 1
      it
    if tr > 0 % check if in refractory period
      v(it) = V_reset; % set voltage to reset
      tr = tr - 1; % reduce running counter of refractory period
    elseif v(it) >= V_th % if voltage over threshold
      rec_spikes(end+1) = it; % record spike event
      v(it) = V_reset; % reset voltage
      tr = tref / dt; % set refractory time
    end
    % Calculate the increment of the membrane potential
    dv = (-g_L * (v(it) - E_L) + Iinj(it)) / tau_m * dt;
    % Update the membrane potential
    v(it + 1) = v(it) + dv;
  end
  % Get spike times in ms
  rec_spikes = rec_spikes * dt;
  rec_v = v;
  rec_sp = rec_spikes;
end