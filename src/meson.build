project(
  'strom',
  'cpp',
  default_options: [
    'cpp_std=c++14',
    'buildtype=release',
    'prefix=/home/basanta/projects/strom/distr',
    'cpp_args=-Wno-deprecated-declarations',
  ],
  version: '1.0',
)

cpp = meson.get_compiler('cpp')
lib_filesystem = cpp.find_library('boost_filesystem', dirs: ['/Users/basanta/lib/static'], required: true)
lib_program_options = cpp.find_library('boost_program_options', dirs: ['/Users/basanta/lib/static'], required: true)
lib_ncl = cpp.find_library('ncl', dirs: ['/Users/basanta/lib/static'], required: true)
lib_beagle = cpp.find_library('hmsbeagle', dirs: ['/usr/local/lib'], required: true)
incl_beagle = include_directories('/usr/local/include/libhmsbeagle-1')
incl_ncl = include_directories('/Users/basanta/include')
incl_eigen = include_directories('/Users/basanta/Documents/libraries/eigen-3.4.0')

executable(
  'strom',
  'main.cpp',
  install: true,
  install_dir: get_option('prefix'),
  dependencies: [lib_beagle, lib_ncl, lib_program_options, lib_filesystem],
  include_directories: [incl_beagle, incl_ncl, incl_eigen],
  link_args: ['-Wl,-rpath,/usr/local/lib'],
)
install_data('strom.conf', install_dir: get_option('prefix'))
install_data('rbcl10.nex', install_dir: get_option('prefix'))
install_data('rbcl10.tre', install_dir: get_option('prefix'))
install_data('rbcl738.nex', install_dir: get_option('prefix'))
install_data('rbcl738nj.tre', install_dir: get_option('prefix'))
install_data('test.tre', install_dir: get_option('prefix'))
install_data('rbcL.nex', install_dir: get_option('prefix'))
install_data('datatest.nex', install_dir: get_option('prefix'))
install_data('rbcLjc.tre', install_dir: get_option('prefix'))
install_data('go.sh', install_dir: '.')
