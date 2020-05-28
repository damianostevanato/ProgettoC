main_iplib: bmp.o ip_lib.o main_iplib.c
	gcc bmp.o ip_lib.o main_iplib.c -o main_iplib -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
bmp.o: bmp.h bmp.c
	gcc -c bmp.c -o bmp.o -Wall -lm
ip_lib.o: ip_lib.h ip_lib.c
	gcc -c ip_lib.c -o ip_lib.o -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
clean: 
	rm bmp.o ip_lib.o main_iplib
	echo Clean Done