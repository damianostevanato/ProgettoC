# ProgettoC
Image processing Project
# Come compilare
# Compilazione normale:
    gcc -c ip_lib.c -o ip_lib.o -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
    gcc ip_lib.o bmp.o main.c -o main -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
# Compilazione per verifica memory leaks:
    gcc -c ip_lib.c -o ip_lib.o -lm -Wall -g -O1
    gcc ip_lib.o bmp.o main.c -o main -lm -Wall -g -O1
    valgrind -v --leak-check=full ./main
