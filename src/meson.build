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
home_dir = '/home/basanta'

boost_dep = dependency('boost', modules: ['filesystem', 'program_options'], required: true)
lib_beagle = cpp.find_library('hmsbeagle', dirs: [home_dir + '/lib/hmsbeagle'], required: true)
lib_ncl = cpp.find_library('ncl', dirs: [home_dir + '/lib/ncl'], required: true)
incl_beagle = include_directories('/home/basanta/include/libhmsbeagle-1')
incl_ncl = include_directories('/home/basanta/include')
incl_eigen = include_directories('/home/basanta/Documents/libraries/eigen-3.4.0')

beagle_lib_dir = home_dir + '/lib/hmsbeagle'
ncl_lib_dir = home_dir + '/lib/ncl'

executable(
  'strom',
  'main.cpp',
  install: true,
  install_dir: get_option('prefix'),
  dependencies: [boost_dep, lib_beagle, lib_ncl],
  include_directories: [incl_beagle, incl_ncl, incl_eigen],
  install_rpath: beagle_lib_dir + ':' + ncl_lib_dir,
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
