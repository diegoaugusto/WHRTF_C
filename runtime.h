#ifndef _RUNTIME_H
#define _RUNTIME_H

#define INIT_VARIABLES	\
	struct timeval inicio;      \
	struct timeval final;       \
	struct timeval resultado;	\
	long int resultado_diff

#define INIT_RUNTIME        \
  gettimeofday(&inicio, NULL)

#define END_RUNTIME           \
  gettimeofday(&final, NULL)

#define PRINT_RUNTIME \
  resultado_diff = (final.tv_usec + 1000000 * final.tv_sec) - (inicio.tv_usec + 1000000 * inicio.tv_sec); \
  resultado.tv_sec = resultado_diff / 1000000; \
  resultado.tv_usec = resultado_diff % 1000000;\
  printf("A execucao durou %lu,%06lu segundos. resultado_diff = %lu\n", resultado.tv_sec, resultado.tv_usec, resultado_diff)

#endif
