project('goldstandard', 'cpp',
        version : '0.0.1',
        license : 'GPL-3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true'])

aahash_dep = declare_dependency(include_directories : 'external/aaHash')

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

src_files = [ 'goldstandard.cpp' ]
deps = [ btllib_dep, threads_dep, openmp_dep, tcmalloc_dep, divsufsort_dep, divsufsort64_dep, sdsl_dep, sequence_dep, z_dep, aahash_dep ]
link_args = [ '-ltcmalloc', '-lsdsl', '-ldivsufsort', '-ldivsufsort64', '-lbtllib', '-lsequence', '-lz' ]

executable('goldstandard',
    src_files,
    dependencies : deps,
    link_args : link_args,
    install : true
)
