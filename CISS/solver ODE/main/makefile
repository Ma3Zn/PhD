USEFUL_FLAGS = 

INC_DIR = ./include/
SRC_DIR = ./source/
OBJ_DIR = ./object/

OBJS := main.o
OBJS := $(addprefix $(OBJ_DIR), $(OBJS))

CC_DEBUG_FLAGS = -DDEBUG
CCFLAGS = -Werror -Wall -Wextra -Wno-unused-parameter -c
#CC = icx
#CC = gcc
CC = icc

.PHONY: build debug clean

build : build_compile

debug : clean debug_compile

build_blis : 
	make build
	make clean
	make build_compile_blis

build_compile : CCFLAGS += -std=c11 -O3
build_compile : $(OBJS)

debug_compile : CCFLAGS += -g3 $(CC_DEBUG_FLAGS)
debug_compile : $(OBJS)

build_compile_blis : CCFLAGS += -O2 -DBLIS_CBLAS
build_compile_blis : $(OBJS)

clean :
	rm $(OBJ_DIR)*.o

$(OBJ_DIR)main.o : $(SRC_DIR)main.c
	$(CC) $(CCFLAGS) $(SRC_DIR)main.c
	mv main.o $@
