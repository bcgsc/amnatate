project('aakomp', 'cpp',
        version : '1.0.0',
        license : 'GPL-3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true'])


threads_dep = dependency('threads')
openmp_dep = dependency('openmp')

compiler = meson.get_compiler('cpp')
tcmalloc_dep = compiler.find_library('tcmalloc')
divsufsort_dep = compiler.find_library('divsufsort')
divsufsort64_dep = compiler.find_library('divsufsort64')
sdsl_dep = compiler.find_library('sdsl')
sequence_dep = compiler.find_library('sequence')
z_dep = compiler.find_library('z')
btllib_dep =  compiler.find_library('btllib')
boost_dep = dependency('boost')

rootpath = meson.project_source_root()

aakomp_src_files = [ 'aakomp.cpp' ]
make_mibf_src_files = [ 'make_mibf.cpp' ]
inc_dirs = [ 'argparse' ]
deps = [ btllib_dep, threads_dep, openmp_dep, tcmalloc_dep, divsufsort_dep, divsufsort64_dep, sdsl_dep, sequence_dep, z_dep, boost_dep ]
link_args = [ '-ltcmalloc', '-lsdsl', '-ldivsufsort', '-ldivsufsort64', '-lbtllib', '-lsequence', '-lz']

executable('aakomp',
    aakomp_src_files,
    dependencies : deps,
    link_args : link_args,
    include_directories : inc_dirs,
    install : true
)

executable('make_mibf',
    make_mibf_src_files,
    dependencies : deps,
    link_args : link_args,
    include_directories : inc_dirs,
    install : true
)

scripts = ['run-aakomp', 'aakomp_score.py']
install_data(scripts, install_dir : 'bin', install_mode : 'rwxr-xr--')
