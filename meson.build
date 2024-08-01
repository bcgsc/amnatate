project('goldstandard', 'cpp',
        version : '0.0.1',
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

rootpath = meson.source_root()

gs_src_files = [ 'goldstandard.cpp' ]
make_mibf_src_files = [ 'make_mibf.cpp' ]
deps = [ btllib_dep, threads_dep, openmp_dep, tcmalloc_dep, divsufsort_dep, divsufsort64_dep, sdsl_dep, sequence_dep, z_dep ]
link_args = [ '-ltcmalloc', '-lsdsl', '-ldivsufsort', '-ldivsufsort64', '-lbtllib', '-lsequence', '-lz' ]

executable('goldstandard',
    gs_src_files,
    dependencies : deps,
    link_args : link_args,
    install : true
)

executable('make_mibf',
    make_mibf_src_files,
    dependencies : deps,
    link_args : link_args,
    install : true
)

scripts = ['goldstandard_run_pipeline.smk']
install_data(scripts, install_dir : 'bin', install_mode : 'rwxr-xr--')