
--- ERROR MESSAGE ---

  File "/home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1130/lib/python3.10/site-packages/mirage/logging/logging_functions.py", line 111, in wrapped
    func(*args, **kwargs)
  File "/home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1130/lib/python3.10/site-packages/mirage/ramp_generator/obs_generator.py", line 1225, in create
    simexp, simzero = self.add_crs_and_noise(self.seed_image, num_integrations=num_integrations)
  File "/home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1130/lib/python3.10/site-packages/mirage/ramp_generator/obs_generator.py", line 222, in add_crs_and_noise
    ramp, rampzero = self.frame_to_ramp(inseed)
  File "/home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1130/lib/python3.10/site-packages/mirage/ramp_generator/obs_generator.py", line 1902, in frame_to_ramp
    poissonsignal = self.do_poisson(deltaframe, self.params['simSignals']['poissonseed'])
  File "/home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1130/lib/python3.10/site-packages/mirage/ramp_generator/obs_generator.py", line 1744, in do_poisson
    newimage = np.random.poisson(signalgain, signalgain.shape).astype(np.float64)
  File "mtrand.pyx", line 3595, in numpy.random.mtrand.RandomState.poisson
  File "_common.pyx", line 877, in numpy.random._common.disc
  File "_common.pyx", line 674, in numpy.random._common.discrete_broadcast_d
  File "_common.pyx", line 408, in numpy.random._common.check_array_constraint
ValueError: lam value too large


