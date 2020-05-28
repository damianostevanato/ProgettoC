# ProgettoC
Image processing Project
# Makefile Instructions
doing "make" command will auto-compile needed files.

doing "make clean" will delete all .o files, recompilation needed after this command is launched.

doing "make bmp.o" / "make ip_lib.o" / "make main_iplib" will compile the specified file ,
because "make" will automatically compile anything that's neccesary even after modifing .c files or .h there's no reason to use this commands.

# Adding new rules to makefile
    target : files needed
        commands
    example:

    we want to compile main_iplib.c
    target = main_iplib

    we need bmp.o and ip_lib.o
    files needed = bmp.o ip_lib.o

    we want to link and compile
    commands: gcc bmp.o ip_lib.o main_iplib.c -o main_iplib 

    Note: bmp.o and ip_lib.o must be created before telling makefile the target "main_iplib"