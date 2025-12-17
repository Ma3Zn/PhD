USEFUL_FLAGS = 

INC_DIR = ./include/
SRC_DIR = ./source/
OBJ_DIR = ./object/

OBJS :=  input.o
OBJS := $(addprefix $(OBJ_DIR), $(OBJS))

CC_DEBUG_FLAGS = -DDEBUG
CCFLAGS = -Werror -Wall -Wextra -Wno-unused-parameter -c
#CC = icx
#CC = gcc
CC = icc

.PHONY: build debug clean

build : build_compile

debug : clean debug_compile

build_compile : CCFLAGS += -std=c11 -O3
build_compile : $(OBJS)

debug_compile : CCFLAGS += -g3 $(CC_DEBUG_FLAGS)
debug_compile : $(OBJS)

clean :
	rm $(OBJ_DIR)*.o

$(OBJ_DIR)input.o : $(SRC_DIR)input.c $(INC_DIR)input.h
	$(CC) $(CCFLAGS) $(SRC_DIR)input.c
	mv input.o $@
