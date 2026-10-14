project('selqie', 'fortran', meson_version: '>=0.50')

f_flags = ['-g', '-Ofast', '-ffree-line-length-none', '-ffixed-line-length-none', '-Wno-maybe-uninitialized']
add_project_arguments(f_flags, language: 'fortran')

fortran_sources = files(
  'constants.f90',
  'determinants.f90',
  'hamiltonian.f90',
  'selection.f90',
  'methods.f90',
  'main.f90',
)

python = import('python').find_installation()
ext_suffix = run_command(
  python, '-c',
  'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX") or ".so")',
  check: true,
).stdout().strip()
f2py_flags = ' '.join(f_flags)
meson.get_compiler('fortran').find_library('lapack')

custom_target(
  'selected_ci',
  input: fortran_sources,
  output: 'selected_ci' + ext_suffix,
  command: [
    python, '-m', 'numpy.f2py', '-c', '-m', 'selected_ci', '@INPUT@',
    'only:', 'solve', ':',
    '--f90flags=' + f2py_flags,
    '--f77flags=' + f2py_flags,
  ] + ['-llapack'],
  build_by_default: true,
)