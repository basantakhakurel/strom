project('strom', 'cpp',
        default_options : ['cpp_std=c++11','cpp_link_args=-static','prefix=/home/basanta/projects/strom/'],
        version : '1.0')

cpp = meson.get_compiler('cpp')
lib_ncl = cpp.find_library('ncl', dirs: ['/home/basanta/lib/static'], required: true)
incl_ncl = include_directories('/home/basanta/include')

executable('strom',
            'main.cpp',
            install: true,
            install_dir: get_option('prefix'),
            dependencies: lib_ncl,
            include_directories: [incl_ncl])

install_data('test.tre', install_dir: get_option('prefix'))
