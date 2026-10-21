project(
  'SELQIE',
  'fortran',
  version : '0.1.0'
)

f_flags = [
  '-g',
  '-Ofast',
  '-ffree-line-length-none',
  '-ffixed-line-length-none',
  '-Wno-maybe-uninitialized'
]

add_project_arguments(
  f_flags,
  language : 'fortran'
)

py = import('python').find_installation(
  pure : false
)

ext_suffix = run_command(
  py,
  '-c',
  'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX") or ".so")',
  check : true
).stdout().strip()

fortran_sources = files(
  'src/hash.f90',
  'src/constants.f90',
  'src/determinants.f90',
  'src/hamiltonian.f90',
  'src/diagonalisers.f90',
  'src/selection.f90',
  'src/methods.f90',
  'src/main.f90'
)

custom_target(
  '_selqie',
  input : fortran_sources,
  output : '_selqie' + ext_suffix,
  command : [
    py,
    '-m',
    'numpy.f2py',
    '-c',
    '-m',
    '_selqie',
    '@INPUT@',
    'only:',
    'solve',
    ':',
    '--f90flags=' + ' '.join(f_flags),
    '--f77flags=' + ' '.join(f_flags),
    '-llapack',
    '-larpack',
    '-lblas'
  ],
  build_by_default : true,
  install : true,
  install_dir : py.get_install_dir() / 'SELQIE'
)

py.install_sources(
  [
    'SELQIE/__init__.py',
    'SELQIE/config.py',
    'SELQIE/utilities.py',
    'SELQIE/selected_ci.py'
  ],
  subdir : 'SELQIE'
)