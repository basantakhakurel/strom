project(
  'strom',
  'cpp',
  default_options: [
    'cpp_std=c++11',
    'prefix=/Users/basantakhakurel/Programming/strom/distr',
  ],
  version: '1.0',
)

cpp = meson.get_compiler('cpp')
lib_filesystem = cpp.find_library(
  'boost_filesystem',
  dirs: ['/Users/basantakhakurel/lib/static'],
  required: true,
)
lib_program_options = cpp.find_library(
  'boost_program_options',
  dirs: ['/Users/basantakhakurel/lib/static'],
  required: true,
)
lib_ncl = cpp.find_library('ncl', dirs: ['/Users/basantakhakurel/lib/static'], required: true)
lib_beagle = cpp.find_library('hmsbeagle', dirs: ['/Users/basantakhakurel/lib'], required: true)
incl_beagle = include_directories('/Users/basantakhakurel/include/libhmsbeagle-1')
incl_ncl = include_directories('/Users/basantakhakurel/include')
incl_boost = include_directories('/Users/basantakhakurel/Documents/libraries/boost_1_86_0/')

executable(
  'strom',
  'main.cpp',
  install: true,
  install_dir: get_option('prefix'),
  dependencies: [lib_beagle, lib_ncl, lib_program_options, lib_filesystem],
  include_directories: [incl_beagle, incl_ncl, incl_boost],
)
install_data('strom.conf', install_dir: get_option('prefix'))
install_data('test.tre', install_dir: get_option('prefix'))
install_data('rbcL.nex', install_dir: get_option('prefix'))
install_data('datatest.nex', install_dir: get_option('prefix'))
install_data('rbcLjc.tre', install_dir: get_option('prefix'))
install_data('go.sh', install_dir: '.')
