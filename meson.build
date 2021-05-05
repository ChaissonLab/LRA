project('lra', 'cpp',
  version : '0.1',
  default_options : ['warning_level=0', 'cpp_std=c++14', 'default_library=static'])

htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])

zlib_dep = dependency('zlib')
pthread_dep = dependency('threads')

executable('lra',
           'lra.cpp',
           install : true,
           dependencies: [htslib_dep, zlib_dep, pthread_dep] )

executable('alchemy2',
           'Alchemy2.cpp',
           install : true,
           dependencies: [htslib_dep, zlib_dep] )

