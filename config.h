#ifndef CONFIG_H
#define CONFIG_H 1

#if __linux__
# define HAVE_FALLOCATE
#else
# undef HAVE_FALLOCATE
#endif

#endif
