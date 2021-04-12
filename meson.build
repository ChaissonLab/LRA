project('lra', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++14'])
htslibdep = dependency('htslib')
zlib_dep = dependency('zlib')
pthread_dep = dependency('threads')

executable('lra',
           'lra.cpp',
           install : true,
           dependencies: [htslibdep, zlib_dep, pthread_dep] )

executable('alchemy2',
           'Alchemy2.cpp',
           install : true,
           dependencies: [htslibdep, zlib_dep] )

